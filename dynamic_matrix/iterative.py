"""Iterative general-formula evaluation with shared, precomputed coefficient sums.

Ternary digit 0: coefficient bit absent / vector bit free.
Digit 1: negative coefficient shift / vector bit fixed to one.
Digit 2: positive coefficient shift / vector bit fixed to zero.
No preparation or evaluation routine in this module recurses.
"""
from dataclasses import dataclass
from functools import reduce
from operator import xor
from typing import Iterable

from .core import Counts, _mul, _power_of_two, _toeplitz_inputs, _values


@dataclass(frozen=True)
class IterativePlan:
    size: int
    coefficients: bytes
    identity_leaves: tuple[int, ...]
    preparation_additions: int

    def matvec(self, vector: Iterable[int], *, counts: Counts | None = None) -> tuple[int, ...]:
        values = _values(vector)
        if len(values) != self.size:
            raise ValueError("vector length must match matrix size")
        work = bytearray(len(self.coefficients))
        work[:self.size] = bytes(values)
        # Expand binary inputs into all distinct subcube sums, in place.
        stride, groups = 1, self.size // 2
        while groups:
            for group in range(groups - 1, -1, -1):
                source, target = group * 2 * stride, group * 3 * stride
                for i in range(stride - 1, -1, -1):
                    a, b = work[source + i], work[source + stride + i]
                    work[target + i] = a ^ b
                    work[target + stride + i] = b
                    work[target + 2 * stride + i] = a
            stride *= 3
            groups //= 2
        # No data-dependent zero/one shortcuts: only structural unit leaves.
        identities = iter(self.identity_leaves)
        next_identity = next(identities, -1)
        for i, coefficient in enumerate(self.coefficients):
            if i == next_identity:
                next_identity = next(identities, -1)
            else:
                work[i] = _mul(coefficient, work[i])
        # Contract ternary products to binary output blocks, in place.
        stride, groups = 1, len(work) // 3
        while groups:
            for group in range(groups):
                source, target = group * 3 * stride, group * 2 * stride
                for i in range(stride):
                    a, b, c = work[source + i], work[source + stride + i], work[source + 2 * stride + i]
                    work[target + i] = a ^ b
                    work[target + stride + i] = a ^ c
            stride *= 2
            groups //= 3
        if counts is not None:
            counts.multiplications += len(self.coefficients) - len(self.identity_leaves)
            counts.additions += 3 * (len(self.coefficients) - self.size)
            counts.identity_copies += len(self.identity_leaves)
        return tuple(work[:self.size])


@dataclass(frozen=True)
class IterativeLayout:
    """Size-only immutable schedule; share across matrices and repeated vectors.

    coefficient_ops is a common-subexpression table: every pair is evaluated once.
    Pairs are canonicalized by operand IDs, never by accidental numeric equality.
    This is a concrete reconstruction of the paper's PARAMS idea, not its lost code.
    """
    size: int
    circulant: bool
    input_count: int
    coefficient_ops: tuple[tuple[int, int], ...]
    coefficient_ids: tuple[int, ...]
    unit_leaves: tuple[int, ...]

    def prepare(self, first_column: Iterable[int], first_row: Iterable[int] | None = None,
                *, unit_row_sum: bool = False) -> IterativePlan:
        if type(unit_row_sum) is not bool:
            raise ValueError("unit_row_sum must be a bool")
        if self.circulant:
            column = _values(first_column)
            if first_row is not None:
                raise ValueError("circulant layouts take a first column only")
            if len(column) != self.size:
                raise ValueError("matrix size must match layout")
            if unit_row_sum and reduce(xor, column) != 1:
                raise ValueError("unit_row_sum=True requires column XOR to equal 1")
            inputs = column
        else:
            if first_row is None or unit_row_sum:
                raise ValueError("Toeplitz layouts require a row and do not accept the unit shortcut")
            column, row = _toeplitz_inputs(first_column, first_row)
            if len(column) != self.size:
                raise ValueError("matrix size must match layout")
            inputs = tuple(reversed(row[1:])) + column
        params = list(inputs) + [0]  # final input slot is a structural zero
        for a, b in self.coefficient_ops:
            params.append(params[a] ^ params[b])
        coefficients = bytes(params[i] for i in self.coefficient_ids)
        identities = self.unit_leaves if unit_row_sum else ()
        if any(coefficients[i] != 1 for i in identities):
            raise AssertionError("unit-row-sum invariant violated")
        return IterativePlan(self.size, coefficients, identities, len(self.coefficient_ops))


def make_iterative_layout(size: int, *, circulant: bool = False) -> IterativeLayout:
    """Build a shared coefficient DAG and ternary leaf order using loops only."""
    if type(size) is not int or type(circulant) is not bool:
        raise ValueError("size must be an integer and circulant must be a bool")
    _power_of_two(size)
    input_count = size if circulant else 2 * size - 1
    zero = input_count
    operations, memo = [], {}

    def intern(a, b):
        if a == b:
            return zero
        if a == zero:
            return b
        if b == zero:
            return a
        pair = (min(a, b), max(a, b))
        if pair not in memo:
            memo[pair] = input_count + 1 + len(operations)
            operations.append(pair)
        return memo[pair]

    current = [d % size for d in range(1-size, size)] if circulant else list(range(input_count))
    width, groups = size, 1
    while width > 1:
        half, child_width = width // 2, width - 1
        following = []
        for group in range(groups):
            start = group * (2 * width - 1)
            a0 = current[start + half:start + half + child_width]
            following.extend(a0)
            following.extend(intern(a0[i], current[start + i]) for i in range(child_width))
            following.extend(intern(a0[i], current[start + width + i]) for i in range(child_width))
        current, width, groups = following, half, groups * 3
    # All nonzero ternary digits mean all coefficient bits are selected.
    unit_leaves = []
    for index in range(len(current)):
        code = index
        for _ in range(size.bit_length() - 1):
            if code % 3 == 0:
                break
            code //= 3
        else:
            unit_leaves.append(index)
    return IterativeLayout(size, circulant, input_count, tuple(operations), tuple(current), tuple(unit_leaves))
