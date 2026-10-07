"""Build a Markdown report from checked-in measurements, without rerunning them."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "benchmarks"


def digest(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def main():
    python_unit = json.loads((RESULTS / "python-unit.json").read_text(encoding="utf-8"))
    lines = [
        "# Local benchmark results", "",
        "Measured on 2026-10-07 (Asia/Taipei). These are new reconstruction measurements, not the manuscript's original results.", "",
        "## Environment and method", "",
        f"- OS: `{python_unit['platform']}`",
        f"- Processor identifier: `{python_unit['processor']}`",
        f"- Python: `{python_unit['python']}`",
        "- C++: Clang 18.1.6 from Zig 0.13.0, native Windows x86-64, C++17, `-O3 -Wall -Wextra -Wpedantic`.",
        "- Python: 7 samples, 3 iterations/sample. C++: 7 samples, 100 iterations/sample.",
        "- Both: seed 20261007, median milliseconds per call, same scalar kernel within each language, counters disabled.",
        "- Fixed matrix/vector per size, warmup, shuffled measurement order; C++ checksum consumption included.",
        "- `apply` reuses a prepared plan; `end-to-end` creates a new plan per product. Matrix preparation and destruction are included in the latter.",
        "- Different language PRNGs yield different benchmark inputs. Shared fixtures are used for cross-language correctness instead.",
        "- Tiny cases are sensitive to clock resolution and overhead. These are single-machine samples, not portable performance guarantees.", "",
        "## N = 128", "",
        "Speedup is baseline time divided by dynamic time. **Below 1 means slower than the baseline.**", "",
        "| Language | Matrix | Baseline ms | Apply ms | End-to-end ms | Apply speedup | End-to-end speedup |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    all_rows = {}
    for language in ("python", "cpp"):
        for kind in ("toeplitz", "circulant", "unit"):
            with (RESULTS / f"{language}-{kind}.csv").open(encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            all_rows[language, kind] = rows
            row = next(r for r in rows if int(r["N"]) == 128)
            lines.append(f"| {language} | {kind} | {float(row['baseline_ms']):.5f} | {float(row['apply_ms']):.5f} | {float(row['end_to_end_ms']):.5f} | {float(row['apply_speedup']):.2f}x | {float(row['end_to_end_speedup']):.2f}x |")
    lines += ["", "In this run, prepared evaluation is faster at N=128 in both languages. C++ end-to-end evaluation is slower at that size because preparation dominates. The C++ evaluator uses one shared scratch buffer; plan construction still allocates a recursion tree. Python's relative speedup does not mean its absolute runtime is lower than C++.",
              "", "## Unit-row-sum circulant scaling", "",
              "| N | Python apply speedup | Python end-to-end speedup | C++ apply speedup | C++ end-to-end speedup |",
              "|---:|---:|---:|---:|---:|"]
    for py, cpp in zip(all_rows["python", "unit"], all_rows["cpp", "unit"]):
        lines.append(f"| {py['N']} | {float(py['apply_speedup']):.2f}x | {float(py['end_to_end_speedup']):.2f}x | {float(cpp['apply_speedup']):.2f}x | {float(cpp['end_to_end_speedup']):.2f}x |")
    lines += ["", "## Artifacts", "",
              "- `python-*.json`: metadata, medians, and all raw samples.",
              "- `python-*.csv`, `cpp-*.csv`: tabular medians and speedups.",
              "- `cpp-*-meta.txt`: compiler and run parameters emitted by the executable. C++ raw samples are not retained by this benchmark runner.",
              "- Regenerate this report with `python tools/summarize_benchmarks.py`.",
              "", "Core source SHA-256 values at report generation:", ""]
    for source in ("dynamic_matrix/core.py", "dynamic_matrix/__main__.py", "cpp/dynamic_matrix.hpp", "cpp/benchmark.cpp"):
        lines.append(f"- `{source}`: `{digest(source)}`")
    (RESULTS / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Wrote benchmarks/README.md")


if __name__ == "__main__":
    main()
