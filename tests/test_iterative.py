import inspect
import random
import unittest
from functools import reduce
from operator import xor

from dynamic_matrix import Counts, gf_mul, make_iterative_layout, prepare_circulant, prepare_toeplitz
from test_core import polynomial_product


def direct_general_formula(column, row, vector, circulant=False):
    """Independent literal SM_A/SM_B oracle, O(N^3); small sizes only.

    For output i and mask j, A=sum a_(p-s), s subset j, p=i&j.
    B=sum b_((~i & j)|t), t subset (~j & (N-1)).
    This writes out the sets; it does not use the implementation's ternary stages.
    """
    size = len(column)
    output = []
    for i in range(size):
        value = 0
        for mask in range(size):
            p = i & mask
            coefficient = 0
            data = 0
            for subset in range(size):
                if subset & mask == subset:
                    offset = p - subset
                    coefficient ^= column[offset % size] if circulant else (column[offset] if offset >= 0 else row[-offset])
                complement = (size - 1) ^ mask
                if subset & complement == subset:
                    data ^= vector[((~i) & mask) | subset]
            value ^= polynomial_product(coefficient, data)
        output.append(value)
    return tuple(output)


class IterativeTests(unittest.TestCase):
    def test_literal_paper_general_formula(self):
        rng = random.Random(992)
        for size in (1, 2, 4, 8, 16):
            column = [rng.randrange(256) for _ in range(size)]
            row = [column[0]] + [rng.randrange(256) for _ in range(size-1)]
            vector = [rng.randrange(256) for _ in range(size)]
            for circular in (False, True):
                layout = make_iterative_layout(size, circulant=circular)
                plan = layout.prepare(column, None if circular else row)
                self.assertEqual(plan.matvec(vector), direct_general_formula(column, row, vector, circular))

    def test_dense_oracle_counts_and_reuse(self):
        rng = random.Random(194)
        for n in range(8):
            size = 2**n
            for circular in (False, True):
                layout = make_iterative_layout(size, circulant=circular)
                self.assertEqual(len(layout.coefficient_ops), len(set(layout.coefficient_ops)))
                for unit in ((False, True) if circular else (False,)):
                    for _ in range(2):
                        column = [rng.randrange(256) for _ in range(size)]
                        if unit:
                            column[-1] = 1 ^ reduce(xor, column[:-1], 0)
                        row = [column[0]] + [rng.randrange(256) for _ in range(size-1)]
                        vector = [rng.randrange(256) for _ in range(size)]
                        plan = layout.prepare(column, None if circular else row, unit_row_sum=unit)
                        recursive = prepare_circulant(column, unit_row_sum=unit) if circular else prepare_toeplitz(column, row)
                        expected = []
                        for i in range(size):
                            expected.append(reduce(xor, (polynomial_product(column[(i-j)%size] if circular else (column[i-j] if i>=j else row[j-i]), vector[j]) for j in range(size)), 0))
                        counts = Counts()
                        self.assertEqual(plan.matvec(vector, counts=counts), tuple(expected))
                        self.assertEqual(plan.matvec(vector), recursive.matvec(vector))
                        self.assertEqual(counts.multiplications, 3**n - (size if unit else 0))
                        self.assertEqual(counts.additions, 3*(3**n-size))
                        self.assertEqual(counts.identity_copies, size if unit else 0)
                        self.assertLessEqual(plan.preparation_additions, recursive.preparation_additions)
                        self.assertEqual(plan.matvec([0]*size), (0,)*size)
                        self.assertEqual(plan.matvec(vector), tuple(expected))

    def test_shared_circulant_coefficients(self):
        layout = make_iterative_layout(128, circulant=True)
        # All N full-row sums share one symbolic coefficient ID.
        self.assertEqual(len({layout.coefficient_ids[i] for i in layout.unit_leaves}), 1)
        self.assertLess(len(layout.coefficient_ops), 3*3**7 - 4*128 + 1)

    def test_invalid_inputs(self):
        for size in (0, 3, True, 1.5):
            with self.assertRaises(ValueError):
                make_iterative_layout(size)
        with self.assertRaises(ValueError):
            make_iterative_layout(2, circulant=True).prepare([1, 1], unit_row_sum=True)
        with self.assertRaises(ValueError):
            make_iterative_layout(2).prepare([1, 2])
        with self.assertRaises(ValueError):
            make_iterative_layout(2).prepare([1, 2], [1, 3], unit_row_sum=True)
        with self.assertRaises(ValueError):
            make_iterative_layout(2, circulant=True).prepare([1])
        with self.assertRaises(ValueError):
            make_iterative_layout(2, circulant=True).prepare([1, 0]).matvec([1])


if __name__ == '__main__':
    unittest.main()
