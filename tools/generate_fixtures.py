"""Regenerate deterministic shared cases from independent polynomial arithmetic."""
import random
import sys
from functools import reduce
from operator import xor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from dynamic_matrix import prepare_circulant, prepare_toeplitz


def product(a, b):
    value = 0
    for bit in range(8):
        if b >> bit & 1:
            value ^= a << bit
    for degree in range(14, 7, -1):
        if value >> degree & 1:
            value ^= 0x11B << (degree - 8)
    return value


def contents():
    rng = random.Random(20261007)
    lines = ["# kind N first_column[N] first_row[N] vector[N] expected[N]; seed=20261007; modulus=0x11B"]
    for exponent in range(8):
        size = 2 ** exponent
        for mode in "TCU":
            for _ in range(3):
                column = [rng.randrange(256) for _ in range(size)]
                if mode == "U":
                    column[-1] = 1 ^ reduce(xor, column[:-1], 0)
                row = [column[0]] + [rng.randrange(256) for _ in range(size-1)] if mode == "T" else column[:1] + column[:0:-1]
                vector = [rng.randrange(256) for _ in range(size)]
                matrix = [[column[i-j] if i >= j else row[j-i] for j in range(size)] for i in range(size)]
                expected = [reduce(xor, (product(a, b) for a, b in zip(r, vector)), 0) for r in matrix]
                plan = prepare_toeplitz(column, row) if mode == "T" else prepare_circulant(column, unit_row_sum=mode == "U")
                if tuple(expected) != plan.matvec(vector):
                    raise AssertionError("fixture verification failed")
                lines.append(" ".join(map(str, [mode, size, *column, *row, *vector, *expected])))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    target = ROOT / "tests" / "fixtures.txt"
    result = contents()
    if "--check" in sys.argv:
        if target.read_text(encoding="utf-8") != result:
            raise SystemExit("fixtures are stale")
        print("72 shared fixtures match Python and independent polynomial oracle")
    else:
        target.write_text(result, encoding="utf-8")
        print("Wrote 72 shared fixtures")
