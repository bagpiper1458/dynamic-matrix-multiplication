# Dynamic Matrix Multiplication

**A reconstruction of Dynamic Multiplication for power-of-two Toeplitz and circulant matrix–vector products over GF(2⁸), in Python and C++17.**

[繁體中文導覽](docs/README.zh-TW.md) · [Derivation and cost model](docs/algorithm.md) · [Reconstruction notes](docs/reconstruction.md)

This project reconstructs the algorithm in the author-provided manuscript *Improve Computing Efficiency of 2ⁿ × 2ⁿ Toeplitz Matrix and Circulant Matrix Multiplication over Galois Fields by Dynamic Multiplication*, by Yang-Bi Jiang, Yan-Wen Chen, Luo-You Jian, and Li-Der Chou. The original source code was lost. This is a new implementation based on the manuscript's block identities, with fresh tests and benchmarks; it is not a recovery of the original code or its experimental environment. No publication venue or status is asserted here.

## What it does

- Computes **matrix–vector** products for N = 2ⁿ, including N = 1.
- Implements the three-recursion decomposition for general Toeplitz matrices and its circulant specialization.
- Supports the unit-row-sum circulant shortcut, with the row sum checked using **XOR**.
- Separates matrix preparation from repeated vector evaluation.
- Counts executed field operations and compares against a conventional O(N²) baseline.
- Includes independent scalar arithmetic checks, shared Python/C++ fixtures, input validation, and reproducible benchmarks.

All elements are bytes representing polynomials modulo **0x11B**, or x⁸ + x⁴ + x³ + x + 1, matching the manuscript's scalar code. Circulant input is the **first column**: `C[i,j] = column[(i-j) % N]`.

| Method | Field multiplications per vector | Field additions per vector |
|---|---:|---:|
| Conventional | N² | N(N−1) |
| Dynamic | 3ⁿ | 3(3ⁿ−2ⁿ) |
| Dynamic, unit-row-sum circulant | 3ⁿ−2ⁿ | 3(3ⁿ−2ⁿ) |

The dynamic addition counts exclude coefficient preparation. This implementation performs **3·3ⁿ − 4·2ⁿ + 1 coefficient XORs** during preparation, plus input checks. These costs and allocation overhead matter for a single product. Fewer field operations do not guarantee a lower runtime. See [the full cost model](docs/algorithm.md).

## Python: no third-party runtime dependencies

Requires Python 3.10 or later. Run from the repository root:

```sh
python -m dynamic_matrix demo
python -m unittest discover -s tests -v
python tools/generate_fixtures.py --check
python -m dynamic_matrix counts --max-exponent 7
```

```python
from dynamic_matrix import Counts, prepare_circulant, prepare_toeplitz

plan = prepare_circulant([2, 1, 1, 3], unit_row_sum=True)
counts = Counts()
result = plan.matvec([0xDB, 0x13, 0x53, 0x45], counts=counts)
assert result == (0x8E, 0x4D, 0xA1, 0xBC)
assert (counts.multiplications, counts.additions) == (5, 15)

# Reuse a plan with a new vector; it owns its prepared coefficients.
assert plan.matvec([1, 0, 0, 0]) == (2, 1, 1, 3)

# T = [[5, 9], [7, 5]]
toeplitz = prepare_toeplitz(first_column=[5, 7], first_row=[5, 9])
print(toeplitz.matvec([11, 13]))
```

`Counts` accumulates when reused. Start a new instance for each measured product. Plans are immutable; input lists are not modified. Inputs must be integer bytes (Python rejects booleans, floats, negatives, and values above 255).

## C++17

The implementation is a standalone header at `cpp/dynamic_matrix.hpp`. The C++ API uses `std::vector<std::uint8_t>`; callers must validate wider numeric input before conversion to bytes.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

With GCC or Clang, CMake is optional:

```sh
c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic cpp/test.cpp -o dm_test
./dm_test tests/fixtures.txt
c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic cpp/benchmark.cpp -o dm_benchmark
./dm_benchmark 7 7 100 circulant
```

On Windows, use `dm_test.exe` / `dm_benchmark.exe`; multi-configuration CMake builds normally place them in `build/Release/`.

```cpp
#include "cpp/dynamic_matrix.hpp"
auto plan = dynamic_matrix::prepare_circulant({2, 1, 1, 3}, true);
auto result = plan.matvec({0xDB, 0x13, 0x53, 0x45});
```

## Benchmarks

See the [recorded local results](benchmarks/README.md). At N=128, this run shows faster prepared evaluation in both languages, while C++ end-to-end evaluation remains slower because matrix preparation dominates.

```sh
python -m dynamic_matrix benchmark --kind circulant --output benchmarks/local-circulant.json
python -m dynamic_matrix benchmark --kind unit-circulant --output benchmarks/local-unit.json
python -m dynamic_matrix benchmark --kind toeplitz --output benchmarks/local-toeplitz.json
```

Python writes JSON with runtime metadata and raw samples, plus CSV. Arguments include `--seed`, `--repeats`, `--iterations`, and `--max-exponent` (0–10). C++ accepts `[max_exponent] [repeats] [iterations] [kind]`, writes CSV to stdout and compiler/measurement metadata to stderr. Redirect these to separate files when saving a run.

Measurements report preparation, prepared application, and preparation-plus-application separately. The conventional baseline operates directly on structured inputs; it does not allocate a dense matrix. Both algorithms use the same eight-step scalar multiplication kernel. No multiplication table, SIMD kernel, or special zero/one shortcut is used in the baseline. The unit-sum dynamic variant only skips the structurally justified identity leaves. Timings include language/runtime overhead and allocations; counters are disabled. Each case is warmed up, measurement order is shuffled, and the median is reported. Input matrices and vectors are fixed within each size, representing repeated workloads. C++ consumes outputs through a volatile checksum to resist dead-code elimination.

Python and C++ random generators differ: compare implementations numerically using the **shared fixtures**, not benchmark seeds. Saved benchmark results apply only to their recorded environment and workload. No original timing figures or speedup percentages are represented as reproduced measurements.

## Tests and provenance

- Exhaustive 256 × 256 scalar multiplication comparison with independent polynomial long division.
- Random Toeplitz/circulant matrices through N = 128, plus zero/identity cases and a known MixColumns vector.
- Exactly 72 checked-in cross-language cases, regenerated with an independent scalar oracle.
- Measured operation counts for general and unit-sum variants, including N = 1.
- Reusable-plan behavior and invalid inputs.
- GitHub Actions configuration for Python and C++ correctness checks; local results and CI results must be distinguished.

The supplied manuscript is not included in this repository. Attribution is preserved above; do not infer a paper publication or license from this reconstruction. A distribution license has not yet been selected. The implementations are research/educational software and are not constant-time cryptographic primitives.
