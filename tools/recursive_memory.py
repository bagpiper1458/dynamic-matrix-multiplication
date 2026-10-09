"""Isolated tracemalloc phases for the recursive Python implementation."""
import gc
import json
import platform
import random
import sys
import tracemalloc
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from dynamic_matrix.core import prepare_circulant, circulant_matvec


def measure(function):
    gc.collect()
    tracemalloc.start()
    value = function()
    gc.collect()  # clear temporary/free-list retention before retained-byte snapshot
    retained, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return value, retained, peak


def main():
    rng = random.Random(20261007)
    rows = []
    for n in range(11):
        size = 2**n
        column = [rng.randrange(256) for _ in range(size)]
        vector = [rng.randrange(256) for _ in range(size)]
        plan, retained, prepare_peak = measure(lambda: prepare_circulant(column))
        actual, result_retained, apply_peak = measure(lambda: plan.matvec(vector))
        one_shot, _, end_peak = measure(lambda: prepare_circulant(column).matvec(vector))
        expected, _, baseline_peak = measure(lambda: circulant_matvec(column, vector))
        if actual != expected or one_shot != expected:
            raise AssertionError("memory run correctness check failed")
        rows.append(dict(n=n, N=size, plan_retained_bytes=retained, prepare_peak_bytes=prepare_peak,
                         apply_peak_extra_bytes=apply_peak, apply_result_retained_bytes=result_retained,
                         end_to_end_peak_bytes=end_peak, baseline_peak_extra_bytes=baseline_peak,
                         node_count=(3**(n+1)-1)//2, recursion_depth=n+1))
        print(f"N={size:4d} retained plan={retained} bytes apply peak extra={apply_peak} bytes", flush=True)
        del plan, actual, one_shot, expected
    output = ROOT / "analysis" / "recursive" / "python-memory.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(dict(python=sys.version, platform=platform.platform(),
        method="tracemalloc; fresh phase per measurement; existing inputs and prepared plan excluded from apply phase; not RSS; tracing disabled in separate timing runs",
        rows=rows), indent=2) + "\n", encoding="utf-8")


if __name__ == '__main__':
    main()
