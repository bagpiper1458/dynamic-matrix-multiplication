import random
import unittest
from functools import reduce
from operator import xor

from dynamic_matrix import (
    Counts, circulant_matvec, dense_matvec, gf_mul,
    prepare_circulant, prepare_toeplitz, toeplitz_matvec,
)


def polynomial_product(a, b):
    # Independent oracle: carryless polynomial product, then long division.
    result = 0
    for shift in range(8):
        if (b >> shift) & 1:
            result ^= a << shift
    for degree in range(14, 7, -1):
        if (result >> degree) & 1:
            result ^= 0x11B << (degree - 8)
    return result


class FieldTests(unittest.TestCase):
    def test_all_65536_scalar_products(self):
        for a in range(256):
            for b in range(256):
                self.assertEqual(gf_mul(a, b), polynomial_product(a, b), (a, b))

    def test_field_identities(self):
        rng = random.Random(61)
        for a in range(256):
            self.assertEqual(gf_mul(a, 0), 0)
            self.assertEqual(gf_mul(a, 1), a)
            if a:
                self.assertTrue(any(gf_mul(a, b) == 1 for b in range(1, 256)))
            b, c = rng.randrange(256), rng.randrange(256)
            self.assertEqual(gf_mul(a, b ^ c), gf_mul(a, b) ^ gf_mul(a, c))

    def test_invalid_elements(self):
        for bad in (-1, 256, 0.5, True, "1"):
            with self.assertRaises(ValueError):
                gf_mul(bad, 1)


class MatrixTests(unittest.TestCase):
    def test_known_mixcolumns_vector_and_orientation(self):
        column = [2, 1, 1, 3]
        vector = [0xDB, 0x13, 0x53, 0x45]
        expected = (0x8E, 0x4D, 0xA1, 0xBC)
        self.assertEqual(circulant_matvec(column, vector), expected)
        for unit in (False, True):
            self.assertEqual(prepare_circulant(column, unit_row_sum=unit).matvec(vector), expected)
        basis = [1, 0, 0, 0]
        self.assertEqual(prepare_circulant(column).matvec(basis), tuple(column))

    def test_random_matrices_against_independent_dense_oracle(self):
        rng = random.Random(20261007)
        for exponent in range(8):
            size = 2 ** exponent
            for _ in range(3):
                column = [rng.randrange(256) for _ in range(size)]
                row = [column[0]] + [rng.randrange(256) for _ in range(size - 1)]
                vector = [rng.randrange(256) for _ in range(size)]
                dense = [[column[i-j] if i >= j else row[j-i] for j in range(size)] for i in range(size)]
                expected = tuple(reduce(xor, (polynomial_product(a, b) for a, b in zip(r, vector)), 0) for r in dense)
                with self.subTest(size=size, kind="toeplitz"):
                    self.assertEqual(dense_matvec(dense, vector), expected)
                    self.assertEqual(toeplitz_matvec(column, row, vector), expected)
                    self.assertEqual(prepare_toeplitz(column, row).matvec(vector), expected)
                for unit in (False, True):
                    if unit:
                        column[-1] = 1 ^ reduce(xor, column[:-1], 0)
                    dense_c = [[column[(i-j) % size] for j in range(size)] for i in range(size)]
                    expected_c = tuple(reduce(xor, (polynomial_product(a, b) for a, b in zip(r, vector)), 0) for r in dense_c)
                    with self.subTest(size=size, kind="circulant", unit=unit):
                        self.assertEqual(circulant_matvec(column, vector), expected_c)
                        self.assertEqual(prepare_circulant(column, unit_row_sum=unit).matvec(vector), expected_c)

    def test_executed_operation_counts(self):
        rng = random.Random(43)
        for exponent in range(8):
            size = 2 ** exponent
            column = [rng.randrange(256) for _ in range(size)]
            column[-1] = 1 ^ reduce(xor, column[:-1], 0)
            row = [column[0]] + [rng.randrange(256) for _ in range(size-1)]
            vector = [rng.randrange(256) for _ in range(size)]
            plans = [(prepare_toeplitz(column, row), 0),
                     (prepare_circulant(column), 0),
                     (prepare_circulant(column, unit_row_sum=True), size)]
            for plan, copies in plans:
                with self.subTest(size=size, copies=copies):
                    counts = Counts()
                    plan.matvec(vector, counts=counts)
                    self.assertEqual(counts.multiplications, 3 ** exponent - copies)
                    self.assertEqual(counts.additions, 3 * (3 ** exponent - size))
                    self.assertEqual(counts.identity_copies, copies)
                    self.assertEqual(plan.preparation_additions, 3 * 3 ** exponent - 4 * size + 1)
            baseline = Counts()
            circulant_matvec(column, vector, counts=baseline)
            self.assertEqual((baseline.multiplications, baseline.additions), (size * size, size * (size-1)))

    def test_zero_identity_and_all_ones(self):
        for size in (1, 2, 4, 16):
            vector = list(range(size))
            for column in ([0] * size, [1] + [0] * (size-1), [1] * size):
                self.assertEqual(prepare_circulant(column).matvec(vector), circulant_matvec(column, vector))
            self.assertEqual(prepare_circulant([1] + [0] * (size-1), unit_row_sum=True).matvec(vector), tuple(vector))

    def test_plan_reuse_and_input_ownership(self):
        column = [2, 1, 1, 3]
        vector = [7, 5, 2, 9]
        plan = prepare_circulant(column, unit_row_sum=True)
        expected = circulant_matvec(column, vector)
        column[0] = 0
        self.assertEqual(plan.matvec(vector), expected)
        self.assertEqual(vector, [7, 5, 2, 9])
        self.assertEqual(plan.matvec([0] * 4), (0,) * 4)
        self.assertEqual(plan.matvec(vector), expected)

    def test_invalid_matrices_and_vectors(self):
        for column in ([], [1, 2, 3], [1, 256], [True]):
            with self.assertRaises(ValueError):
                prepare_circulant(column)
        for column, row in (([1, 2], [3, 4]), ([1, 2], [1]), ([1], [])):
            with self.assertRaises(ValueError):
                prepare_toeplitz(column, row)
        with self.assertRaises(ValueError):
            prepare_circulant([1, 1], unit_row_sum=True)
        with self.assertRaises(ValueError):
            prepare_circulant([1], unit_row_sum="yes")
        with self.assertRaises(ValueError):
            prepare_circulant([1, 0]).matvec([1])
        with self.assertRaises(ValueError):
            dense_matvec([[1], [1, 2]], [1])


if __name__ == "__main__":
    unittest.main()
