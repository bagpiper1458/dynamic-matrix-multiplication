"""Reproducible demonstrations and benchmarks: python -m dynamic_matrix --help."""

import argparse
import csv
import json
import platform
import random
import statistics
import sys
from datetime import datetime, timezone
from functools import reduce
from operator import xor
from pathlib import Path
from time import perf_counter_ns

from .core import Counts, circulant_matvec, prepare_circulant, prepare_toeplitz, toeplitz_matvec


def demo():
    column, vector = [2, 1, 1, 3], [0xDB, 0x13, 0x53, 0x45]
    for unit in (False, True):
        counts = Counts()
        plan = prepare_circulant(column, unit_row_sum=unit)
        result = plan.matvec(vector, counts=counts)
        print(f"unit_row_sum={unit}: {[f'{v:02x}' for v in result]}, {counts}")


def operation_counts(max_exponent):
    print("n,N,baseline_M,baseline_A,dynamic_M,dynamic_A,unit_M,preparation_A")
    for n in range(max_exponent + 1):
        size = 2 ** n
        column = [1] + [0] * (size - 1)
        normal, unit = Counts(), Counts()
        plan = prepare_circulant(column)
        plan.matvec(column, counts=normal)
        prepare_circulant(column, unit_row_sum=True).matvec(column, counts=unit)
        print(f"{n},{size},{size*size},{size*(size-1)},{normal.multiplications},"
              f"{normal.additions},{unit.multiplications},{plan.preparation_additions}")


def _sample_ms(function, iterations):
    start = perf_counter_ns()
    for _ in range(iterations):
        function()
    return (perf_counter_ns() - start) / iterations / 1e6


def benchmark(args):
    rng = random.Random(args.seed)
    rows = []
    for n in range(args.max_exponent + 1):
        size = 2 ** n
        column = [rng.randrange(256) for _ in range(size)]
        vector = [rng.randrange(256) for _ in range(size)]
        if args.kind == "unit-circulant":
            column[-1] = 1 ^ reduce(xor, column[:-1], 0)
        if args.kind == "toeplitz":
            row = [column[0]] + [rng.randrange(256) for _ in range(size - 1)]
            prepare = lambda: prepare_toeplitz(column, row)
            conventional = lambda: toeplitz_matvec(column, row, vector)
        else:
            prepare = lambda: prepare_circulant(column, unit_row_sum=args.kind == "unit-circulant")
            conventional = lambda: circulant_matvec(column, vector)
        plan = prepare()
        expected = conventional()
        if plan.matvec(vector) != expected:
            raise AssertionError("benchmark correctness check failed")
        counts = Counts()
        plan.matvec(vector, counts=counts)
        cases = {
            "baseline": conventional,
            "prepare": prepare,
            "apply": lambda: plan.matvec(vector),
            "end_to_end": lambda: prepare().matvec(vector),
        }
        samples = {name: [] for name in cases}
        for function in cases.values():
            function()  # untimed warmup
        for _ in range(args.repeats):
            order = list(cases)
            rng.shuffle(order)
            for name in order:
                samples[name].append(_sample_ms(cases[name], args.iterations))
        medians = {name: statistics.median(values) for name, values in samples.items()}
        record = dict(n=n, N=size, kind=args.kind, multiplications=counts.multiplications,
                      additions=counts.additions, preparation_additions=plan.preparation_additions,
                      **{name + "_ms": value for name, value in medians.items()},
                      apply_speedup=medians["baseline"] / medians["apply"],
                      end_to_end_speedup=medians["baseline"] / medians["end_to_end"],
                      samples_ms=samples)
        rows.append(record)
        print(f"N={size:4d} baseline={medians['baseline']:.4f}ms "
              f"apply={medians['apply']:.4f}ms end-to-end={medians['end_to_end']:.4f}ms "
              f"speedup={record['apply_speedup']:.2f}x/{record['end_to_end_speedup']:.2f}x")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = dict(
        generated_utc=datetime.now(timezone.utc).isoformat(),
        python=sys.version, platform=platform.platform(), processor=platform.processor(),
        seed=args.seed, repeats=args.repeats, iterations=args.iterations,
        field_polynomial="0x11B", implementation="Python reference reconstruction",
        timing="median ms per call; counters disabled; randomized case order; same scalar kernel",
        rows=rows,
    )
    output.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    with output.with_suffix(".csv").open("w", newline="", encoding="utf-8") as handle:
        fields = [key for key in rows[0] if key != "samples_ms"]
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved {output} and {output.with_suffix('.csv')}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("demo")
    counts = sub.add_parser("counts")
    counts.add_argument("--max-exponent", type=int, default=7, choices=range(11))
    bench = sub.add_parser("benchmark")
    bench.add_argument("--kind", choices=("toeplitz", "circulant", "unit-circulant"), default="circulant")
    bench.add_argument("--max-exponent", type=int, default=7, choices=range(11))
    bench.add_argument("--repeats", type=int, default=7)
    bench.add_argument("--iterations", type=int, default=3)
    bench.add_argument("--seed", type=int, default=20261007)
    bench.add_argument("--output", default="benchmarks/local.json")
    args = parser.parse_args()
    if args.command == "demo":
        demo()
    elif args.command == "counts":
        operation_counts(args.max_exponent)
    else:
        if args.repeats < 1 or args.iterations < 1:
            parser.error("repeats and iterations must be positive")
        benchmark(args)


if __name__ == "__main__":
    main()
