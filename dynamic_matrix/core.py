"""Paper-derived block recursion; coefficients use polynomial basis modulo 0x11B.

For a Toeplitz block matrix [[A0, A1], [A2, A0]], compute
F=A0(x0 XOR x1), G=(A0 XOR A1)x1, H=(A0 XOR A2)x0.
The output is [F XOR G, F XOR H]. See docs/algorithm.md.
"""

from dataclasses import dataclass
from functools import reduce
from operator import xor
from typing import Iterable


def _byte(value: int) -> None:
    if type(value) is not int or not 0 <= value <= 255:
        raise ValueError("field elements must be integers in [0, 255]")


def _values(values: Iterable[int]) -> tuple[int, ...]:
    result = tuple(values)
    for value in result:
        _byte(value)
    return result


def _power_of_two(size: int) -> None:
    if size < 1 or size & (size - 1):
        raise ValueError("matrix size must be a positive power of two")


def _mul(a: int, b: int) -> int:
    """Internal scalar kernel shared by both algorithms; inputs already checked."""
    product = 0
    for _ in range(8):
        if b & 1:
            product ^= a
        a = ((a << 1) ^ (0x11B if a & 0x80 else 0)) & 0xFF
        b >>= 1
    return product


def gf_mul(a: int, b: int) -> int:
    """Multiply byte polynomials in GF(2^8), modulo x^8+x^4+x^3+x+1."""
    _byte(a)
    _byte(b)
    return _mul(a, b)


@dataclass
class Counts:
    """Executed field operations; scalar-kernel bit operations are excluded.

    A supplied counter accumulates across calls. Use a fresh instance per product.
    Copies count only the identity leaves enabled by the unit-row-sum option.
    """

    multiplications: int = 0
    additions: int = 0
    identity_copies: int = 0


@dataclass(frozen=True)
class _Node:
    size: int
    coefficient: int = 0
    identity: bool = False
    children: tuple["_Node", ...] = ()


def _compile(diagonals: tuple[int, ...], size: int, unit: bool) -> tuple[_Node, int]:
    # diagonals[size-1+d] is the coefficient on row-minus-column diagonal d.
    if size == 1:
        if unit and diagonals[0] != 1:
            raise AssertionError("unit-row-sum invariant violated")
        return _Node(1, diagonals[0], unit), 0
    half = size // 2
    middle = size - 1
    offsets = range(1 - half, half)
    a0 = tuple(diagonals[middle + d] for d in offsets)
    a1 = tuple(diagonals[middle + d - half] for d in offsets)
    a2 = tuple(diagonals[middle + d + half] for d in offsets)
    g_coeff = tuple(a ^ b for a, b in zip(a0, a1))
    h_coeff = tuple(a ^ b for a, b in zip(a0, a2))
    # If the parent is unit-sum circulant, G and H are also unit-sum circulant.
    # F is merely Toeplitz: never propagate the shortcut into F.
    f, f_cost = _compile(a0, half, False)
    g, g_cost = _compile(g_coeff, half, unit)
    h, h_cost = _compile(h_coeff, half, unit)
    return _Node(size, children=(f, g, h)), 2 * (size - 1) + f_cost + g_cost + h_cost


def _evaluate(node: _Node, vector: tuple[int, ...], counts: Counts | None) -> tuple[int, ...]:
    if node.size == 1:
        if node.identity:
            if counts is not None:
                counts.identity_copies += 1
            return vector
        if counts is not None:
            counts.multiplications += 1
        return (_mul(node.coefficient, vector[0]),)
    half = node.size // 2
    left, right = vector[:half], vector[half:]
    mixed = tuple(a ^ b for a, b in zip(left, right))
    f = _evaluate(node.children[0], mixed, counts)
    g = _evaluate(node.children[1], right, counts)
    h = _evaluate(node.children[2], left, counts)
    if counts is not None:
        counts.additions += 3 * half  # mixed, upper result, lower result
    return tuple(a ^ b for a, b in zip(f, g)) + tuple(a ^ b for a, b in zip(f, h))


@dataclass(frozen=True)
class Plan:
    """Reusable immutable matrix plan. Create through prepare_* functions.

    preparation_additions counts coefficient XORs in this implementation.
    It excludes input validation (including the row-sum check) and allocations.
    """

    size: int
    preparation_additions: int
    _root: _Node

    def matvec(self, vector: Iterable[int], *, counts: Counts | None = None) -> tuple[int, ...]:
        values = _values(vector)
        if len(values) != self.size:
            raise ValueError("vector length must match matrix size")
        return _evaluate(self._root, values, counts)


def _toeplitz_inputs(column: Iterable[int], row: Iterable[int]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    column, row = _values(column), _values(row)
    _power_of_two(len(column))
    if len(row) != len(column) or row[0] != column[0]:
        raise ValueError("first row and column must have equal lengths and equal first elements")
    return column, row


def prepare_toeplitz(first_column: Iterable[int], first_row: Iterable[int]) -> Plan:
    """Precompute a power-of-two Toeplitz matrix from its first column and row."""
    column, row = _toeplitz_inputs(first_column, first_row)
    root, additions = _compile(tuple(reversed(row[1:])) + column, len(column), False)
    return Plan(len(column), additions, root)


def prepare_circulant(first_column: Iterable[int], *, unit_row_sum: bool = False) -> Plan:
    """Prepare C[i,j] = first_column[(i-j) % N]. Row sums are field XORs.

    With unit_row_sum=True, require XOR(first_column)==1 and skip exactly N
    structurally identified multiplications by one, including N=1.
    """
    column = _values(first_column)
    _power_of_two(len(column))
    if type(unit_row_sum) is not bool:
        raise ValueError("unit_row_sum must be a bool")
    if unit_row_sum and reduce(xor, column) != 1:
        raise ValueError("unit_row_sum=True requires XOR of the first column to equal 1")
    row = column[:1] + tuple(reversed(column[1:]))
    root, additions = _compile(tuple(reversed(row[1:])) + column, len(column), unit_row_sum)
    return Plan(len(column), additions, root)


def _dot(coefficients: Iterable[int], vector: tuple[int, ...], counts: Counts | None) -> int:
    products = iter(zip(coefficients, vector))
    a, b = next(products)
    total = _mul(a, b)
    for a, b in products:
        total ^= _mul(a, b)
    if counts is not None:
        counts.multiplications += len(vector)
        counts.additions += len(vector) - 1
    return total


def toeplitz_matvec(first_column: Iterable[int], first_row: Iterable[int], vector: Iterable[int], *, counts: Counts | None = None) -> tuple[int, ...]:
    """Conventional O(N^2) baseline using structured input, without a dense allocation."""
    column, row = _toeplitz_inputs(first_column, first_row)
    values = _values(vector)
    if len(values) != len(column):
        raise ValueError("vector length must match matrix size")
    return tuple(_dot((column[i-j] if i >= j else row[j-i] for j in range(len(values))), values, counts)
                 for i in range(len(values)))


def circulant_matvec(first_column: Iterable[int], vector: Iterable[int], *, counts: Counts | None = None) -> tuple[int, ...]:
    """Conventional O(N^2) circulant baseline, C[i,j] = column[(i-j) % N]."""
    column, values = _values(first_column), _values(vector)
    _power_of_two(len(column))
    if len(values) != len(column):
        raise ValueError("vector length must match matrix size")
    size = len(column)
    return tuple(_dot((column[(i-j) % size] for j in range(size)), values, counts) for i in range(size))


def dense_matvec(matrix: Iterable[Iterable[int]], vector: Iterable[int], *, counts: Counts | None = None) -> tuple[int, ...]:
    """Independent dense indexing baseline, also accepting non-power-of-two sizes."""
    rows, values = tuple(_values(row) for row in matrix), _values(vector)
    if not rows or not values or any(len(row) != len(values) for row in rows):
        raise ValueError("matrix must be nonempty and all rows must match vector length")
    return tuple(_dot(row, values, counts) for row in rows)
