"""Small-size symbolic/numeric audit of a concrete PARAMS interpretation.

Uses the author's recovered lowbit set generator and aligned square blocks.
This is an explanatory reference, not a production implementation or benchmark.
"""
import json
import random
import sys
from functools import reduce
from operator import xor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from dynamic_matrix import make_iterative_layout


def polynomial_product(a, b):
    product = 0
    for bit in range(8):
        if (b >> bit) & 1:
            product ^= a << bit
    for degree in range(14, 7, -1):
        if (product >> degree) & 1:
            product ^= 0x11B << (degree-8)
    return product


def recovered_sets(size):
    """Old C/L/CC parent/shift, with a natural N=1 base case added."""
    rows = []
    for i in range(size):
        a, b = [frozenset({0})], [frozenset({size-1-i})]
        for j in range(1, size):
            step = j & -j
            shift = step if i & step else -step
            a.append(a[j-step] | {index+shift for index in a[j-step]})
            b.append(b[j-step] | {index+shift for index in b[j-step]})
        rows.append(list(zip(a, reversed(b))))
    return rows


def pos(level, row_start, column_start):
    mask = column_start >> level
    return level, mask, (row_start >> level) & mask


def ternary_id(mask, selected, bits):
    return sum((0 if not (mask >> bit) & 1 else 2 if (selected >> bit) & 1 else 1)*3**bit
               for bit in range(bits))


def audit(exponent, kind):
    size, leaves = 2**exponent, 3**exponent
    circular, unit = kind != "toeplitz", kind == "unit-circulant"
    formulas = recovered_sets(size)
    rng = random.Random(20261009+exponent)
    column = [rng.randrange(256) for _ in range(size)]
    row = [column[0]] + [rng.randrange(256) for _ in range(size-1)]
    vector = [rng.randrange(256) for _ in range(size)]
    if unit:
        column[-1] = 1 ^ reduce(xor, column[:-1], 0)
        row[0] = column[0]

    products, vectors, coefficient_forms = {}, {}, set()
    for i, terms in enumerate(formulas):
        for mask, (a, b) in enumerate(terms):
            key = pos(0, i, mask)
            a = frozenset(index % size if circular else index for index in a)
            selected = i & mask
            free, fixed = (size-1)^mask, mask^selected
            assert b == frozenset(fixed | s for s in range(size) if s & free == s)
            if key in products:
                assert products[key] == (a, b)
            products[key] = a, b
            vectors[free, fixed] = b
            coefficient_forms.add(a)
    assert len(products) == len(set(products.values())) == len(vectors) == leaves

    # Every non-singleton subcube sum is one XOR of two smaller subcubes.
    vector_values, vector_xors = {}, 0
    for free, fixed in sorted(vectors, key=lambda k: k[0].bit_count()):
        if free == 0:
            vector_values[free, fixed] = vector[fixed]
        else:
            bit = 1 << (free.bit_length()-1)
            left, right = (free^bit, fixed), (free^bit, fixed|bit)
            assert vectors[free, fixed] == vectors[left] | vectors[right]
            vector_values[free, fixed] = vector_values[left] ^ vector_values[right]
            vector_xors += 1
    assert vector_xors == leaves-size

    def coefficient(index):
        return column[index % size] if circular else column[index] if index >= 0 else row[-index]

    block_values, block_symbols = {}, {}
    multiplies = copies = 0
    for key, (a, b) in products.items():
        _, mask, selected = key
        av = reduce(xor, map(coefficient, a), 0)
        bv = vector_values[(size-1)^mask, mask^selected]
        if unit and mask == size-1:
            assert av == 1
            value = bv
            copies += 1
        else:
            value = polynomial_product(av, bv)
            multiplies += 1
        block_values[key] = (value,)
        # Keep unrestricted coefficient symbols, even for the unit numeric case.
        block_symbols[key] = (frozenset((ai, bj) for ai in a for bj in b),)

    levels = [{"level": 0, "block_side": 1, "unique_blocks": leaves,
               "value_bytes_if_packed": leaves, "new_output_xors": 0}]
    block_xors = 0
    for level in range(1, exponent+1):
        side, half = 2**level, 2**(level-1)
        signatures = {}
        for r in range(0, size, side):
            for c in range(0, size, side):
                key = pos(level, r, c)
                signature = tuple(pos(0, i, j) for i in range(r, r+side) for j in range(c, c+side))
                if key in signatures:
                    assert signatures[key] == signature
                    continue
                signatures[key] = signature
                f = pos(level-1, r, c)
                g = pos(level-1, r, c+half)
                h = pos(level-1, r+half, c+half)
                assert f == pos(level-1, r+half, c)
                # Explicit closed-form child keys, independent of grid positions.
                _, mask, selected = key
                assert (f, g, h) == ((level-1, 2*mask, 2*selected),
                                     (level-1, 2*mask+1, 2*selected),
                                     (level-1, 2*mask+1, 2*selected+1))
                tid = ternary_id(mask, selected, exponent-level)
                assert [ternary_id(k[1], k[2], exponent-level+1) for k in (f,g,h)] == [3*tid,3*tid+1,3*tid+2]
                block_values[key] = tuple(x^y for x,y in zip(block_values[f], block_values[g])) + \
                                    tuple(x^y for x,y in zip(block_values[f], block_values[h]))
                block_symbols[key] = tuple(x^y for x,y in zip(block_symbols[f], block_symbols[g])) + \
                                     tuple(x^y for x,y in zip(block_symbols[f], block_symbols[h]))
                block_xors += side
        count = len(signatures)
        assert count == 3**(exponent-level)
        levels.append({"level": level, "block_side": side, "unique_blocks": count,
                       "value_bytes_if_packed": side*count, "new_output_xors": side*count})

    expected = tuple(reduce(xor, (polynomial_product(coefficient(i-j), vector[j])
                                 for j in range(size)), 0) for i in range(size))
    root = exponent, 0, 0
    assert block_values[root] == expected
    assert block_symbols[root] == tuple(frozenset(((i-j)%size if circular else i-j, j)
                                                  for j in range(size)) for i in range(size))
    plan = make_iterative_layout(size, circulant=circular).prepare(
        column, None if circular else row, unit_row_sum=unit)
    assert plan.matvec(vector) == expected
    assert multiplies == leaves-(size if unit else 0)
    assert copies == (size if unit else 0)
    assert block_xors == 2*(leaves-size)
    assert sum(k["value_bytes_if_packed"] for k in levels) == 3*leaves-2*size
    result = {"exponent": exponent, "N": size, "kind": kind,
              "formula_product_occurrences": size*size, "unique_products": leaves,
              "unique_vector_forms": leaves, "unique_coefficient_forms": len(coefficient_forms),
              "vector_xors": vector_xors, "output_xors": block_xors,
              "field_multiplies": multiplies, "identity_copies": copies,
              "levels": levels, "symbolic_identity": True, "numeric_identity": True,
              "matches_existing_iterative": True}
    if size == 4:
        result["product_id_grid"] = [[ternary_id(j, i&j, exponent) for j in range(size)] for i in range(size)]
    return result


def main():
    results = [audit(p, kind) for p in range(7) for kind in ("toeplitz", "circulant", "unit-circulant")]
    destination = ROOT / "analysis/params-derivation"
    destination.mkdir(parents=True, exist_ok=True)
    (destination / "audit.json").write_text(json.dumps({
        "description": "Finite-size symbolic and numeric checks of derived PARAMS semantics; not a speed/memory benchmark.",
        "sizes": [2**p for p in range(7)], "cases": results,
    }, indent=2), encoding="utf-8")
    for case in results:
        if case["kind"] == "toeplitz":
            print(f"N={case['N']}: {case['unique_products']} products, {case['vector_xors']} vector XORs, "
                  f"{case['output_xors']} output XORs; block keys/symbolic identity/numeric identity PASS")
    print(f"PASS: {len(results)} size/kind cases; all levels and three-child block relations checked")


if __name__ == "__main__":
    main()
