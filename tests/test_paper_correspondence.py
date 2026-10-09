"""Independent Algorithm 1/2 set generation, not the ternary implementation.

Conventions follow the manuscript's displayed 4x4 expansion: SM_B begins at
N-1-i and its columns are reversed when pairing with SM_A. Manuscript prose
also contains b_i; this interpretation is explicit rather than silently chosen.
"""
import random
import unittest
from functools import reduce
from operator import xor

from dynamic_matrix import make_iterative_layout, prepare_circulant, prepare_toeplitz
from test_core import polynomial_product


def paper_sets(size):
    result = []
    for i in range(size):
        a, b = [frozenset({0})], [frozenset({size-1-i})]
        for j in range(1, size):
            step = 1 << (j.bit_length()-1)
            shift = step if i & step else -step
            a.append(a[j-step] | {index+shift for index in a[j-step]})
            b.append(b[j-step] | {index+shift for index in b[j-step]})
        result.append(list(zip(a, reversed(b))))
    return result


class PaperCorrespondenceTests(unittest.TestCase):
    def test_symbolic_coefficient_mapping_and_bilinear_identity(self):
        # These certificates cover every matrix/vector byte assignment at each
        # tested dimension, because scalar bilinearity expands into a[k]*b[j].
        for size in (1, 2, 4, 8, 16, 32):
            formulas = paper_sets(size)
            for circular in (False, True):
                layout = make_iterative_layout(size, circulant=circular)
                inputs = range(size) if circular else range(1-size, size)
                forms = [frozenset({i}) for i in inputs] + [frozenset()]
                for a, b in layout.coefficient_ops:
                    forms.append(forms[a] ^ forms[b])
                products_seen = set()
                for i, terms in enumerate(formulas):
                    expanded = set()
                    for mask, (a, b) in enumerate(terms):
                        a = frozenset(offset % size if circular else offset for offset in a)
                        leaf = sum((0 if not (mask >> bit) & 1 else 2 if (i >> bit) & 1 else 1)*3**bit
                                   for bit in range(size.bit_length()-1))
                        self.assertEqual(forms[layout.coefficient_ids[leaf]], a)
                        self.assertTrue(all(0 <= index < size for index in b))
                        products_seen.add((a, b))
                        expanded.symmetric_difference_update((offset, index) for offset in a for index in b)
                    expected = {((i-j) % size if circular else i-j, j) for j in range(size)}
                    self.assertEqual(expanded, expected)
                self.assertEqual(len(products_seen), 3**(size.bit_length()-1))

    def test_numeric_literal_algorithms_against_both_python_paths(self):
        rng = random.Random(20261009)
        for size in (1, 2, 4, 8, 16):
            formulas = paper_sets(size)
            for circular, unit in ((False, False), (True, False), (True, True)):
                column = [rng.randrange(256) for _ in range(size)]
                if unit:
                    column[-1] = 1 ^ reduce(xor, column[:-1], 0)
                row = [column[0]] + [rng.randrange(256) for _ in range(size-1)]
                vector = [rng.randrange(256) for _ in range(size)]
                def coefficient(offset):
                    return column[offset % size] if circular else column[offset] if offset >= 0 else row[-offset]
                expected = tuple(reduce(xor, (
                    polynomial_product(reduce(xor, map(coefficient, a), 0), reduce(xor, (vector[j] for j in b), 0))
                    for a, b in terms), 0) for terms in formulas)
                recursive = prepare_circulant(column, unit_row_sum=unit) if circular else prepare_toeplitz(column, row)
                iterative = make_iterative_layout(size, circulant=circular).prepare(
                    column, None if circular else row, unit_row_sum=unit)
                self.assertEqual(recursive.matvec(vector), expected)
                self.assertEqual(iterative.matvec(vector), expected)


if __name__ == "__main__":
    unittest.main()
