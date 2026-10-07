# Local benchmark results

Measured on 2026-10-07 (Asia/Taipei). These are new reconstruction measurements, not the manuscript's original results.

## Environment and method

- OS: `Windows-11-10.0.26200-SP0`
- Processor identifier: `AMD64 Family 25 Model 117 Stepping 2, AuthenticAMD`
- Python: `3.12.14 (main, Aug 25 2026, 14:01:42) [MSC v.1944 64 bit (AMD64)]`
- C++: Clang 18.1.6 from Zig 0.13.0, native Windows x86-64, C++17, `-O3 -Wall -Wextra -Wpedantic`.
- Python: 7 samples, 3 iterations/sample. C++: 7 samples, 100 iterations/sample.
- Both: seed 20261007, median milliseconds per call, same scalar kernel within each language, counters disabled.
- Fixed matrix/vector per size, warmup, shuffled measurement order; C++ checksum consumption included.
- `apply` reuses a prepared plan; `end-to-end` creates a new plan per product. Matrix preparation and destruction are included in the latter.
- Different language PRNGs yield different benchmark inputs. Shared fixtures are used for cross-language correctness instead.
- Tiny cases are sensitive to clock resolution and overhead. These are single-machine samples, not portable performance guarantees.

## N = 128

Speedup is baseline time divided by dynamic time. **Below 1 means slower than the baseline.**

| Language | Matrix | Baseline ms | Apply ms | End-to-end ms | Apply speedup | End-to-end speedup |
|---|---|---:|---:|---:|---:|---:|
| python | toeplitz | 13.77200 | 3.24133 | 7.02463 | 4.25x | 1.96x |
| python | circulant | 14.40787 | 3.40920 | 7.20213 | 4.23x | 2.00x |
| python | unit | 14.48653 | 3.22297 | 7.00987 | 4.49x | 2.07x |
| cpp | toeplitz | 0.05614 | 0.01404 | 0.17590 | 4.00x | 0.32x |
| cpp | circulant | 0.05921 | 0.01407 | 0.17720 | 4.21x | 0.33x |
| cpp | unit | 0.05930 | 0.01397 | 0.17939 | 4.25x | 0.33x |

In this run, prepared evaluation is faster at N=128 in both languages. C++ end-to-end evaluation is slower at that size because preparation dominates. The C++ evaluator uses one shared scratch buffer; plan construction still allocates a recursion tree. Python's relative speedup does not mean its absolute runtime is lower than C++.

## Unit-row-sum circulant scaling

| N | Python apply speedup | Python end-to-end speedup | C++ apply speedup | C++ end-to-end speedup |
|---:|---:|---:|---:|---:|
| 1 | 6.10x | 1.11x | 1.12x | 0.18x |
| 2 | 1.94x | 0.68x | 0.70x | 0.14x |
| 4 | 1.74x | 0.69x | 0.95x | 0.11x |
| 8 | 1.78x | 0.76x | 1.37x | 0.11x |
| 16 | 2.10x | 0.94x | 1.92x | 0.14x |
| 32 | 2.67x | 1.21x | 2.57x | 0.19x |
| 64 | 3.37x | 1.53x | 3.30x | 0.25x |
| 128 | 4.49x | 2.07x | 4.25x | 0.33x |

## Artifacts

- `python-*.json`: metadata, medians, and all raw samples.
- `python-*.csv`, `cpp-*.csv`: tabular medians and speedups.
- `cpp-*-meta.txt`: compiler and run parameters emitted by the executable. C++ raw samples are not retained by this benchmark runner.
- Regenerate this report with `python tools/summarize_benchmarks.py`.

Core source SHA-256 values at report generation:

- `dynamic_matrix/core.py`: `9c6534b33e7b64ecbe60448d15294ee9221fc5424371b4e61046bdf8cc0022d3`
- `dynamic_matrix/__main__.py`: `7b1a0c09e9bf7987808fdd7063b5ff0ff52993499b7e9cdc053969521b3b8faf`
- `cpp/dynamic_matrix.hpp`: `84744f3b392ab1305293e7e9b2f847466156b0b5950290586281384e6f7daab8`
- `cpp/benchmark.cpp`: `5cb0489dcaa382184c24251a1312a3c6f7210afb39a6ea976aaf9e3f2b6d115d`
