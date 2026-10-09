# Compact C++ implementations and fair comparison

The earlier code passed correctness checks, but its recursive plan retained a pointer tree while its iterative plan used a compact coefficient array. That comparison mixed traversal costs with storage choices. `cpp/packed_matrix.hpp` now provides improved brute force and two evaluators backed by the **same coefficient representation**, while preserving the earlier implementations as references.

This is an engineering optimization under the paper's decomposition, not a proof of global optimality. Python remains the independent formula/reference implementation; the new performance variants are C++17.

## What changed

- **Brute force:** circulant uses `(i-j) & (N-1)` instead of integer remainder, valid because sizes are checked to be powers of two. Toeplitz splits each row's loop at the diagonal, avoiding a coefficient-selection branch per element. Both still execute N² field multiplications.
- **Recursive preparation:** recursively fills a flat leaf-coefficient array. Child coefficient diagonals use a preallocated scratch area; there are no per-node heap allocations or pointer tree.
- **Iterative preparation:** starts from the 2N−1 diagonals in the coefficient array, then loops through F/G/H levels. It processes source groups backwards, saving one source diagonal group before overwriting it. No recursion, symbolic hash table, or retained layout is required.
- **Recursive application:** recursively traverses the flat coefficient array with known child offsets; a fused 2×2 base case avoids three scalar recursive calls. It reuses scratch across child calls.
- **Iterative application:** preserves the proven shared-vector expansion and PARAMS block contraction. General multiplication has no per-element identity branch. For unit-sum circulants, a binary counter enumerates the all-nonzero ternary leaf IDs and multiplies the intervening ranges; the N structural identities are copied without field multiplication.
- **Inner loops:** iterative preparation, expansion, multiplication and contraction acquire the vector data pointers before their inner loops, making contiguous access explicit to the compiler.
- **Allocations:** every method offers both an allocating wrapper and a caller-owned output/workspace API. The latter performs zero heap allocations per application; its retained buffers are still counted in memory reports.

All paths use the same eight-step `gf_mul` modulo 0x11B. There are no multiplication lookup tables or explicit SIMD kernels. Compiler auto-vectorization is allowed for all implementations. None of these changes changes the definition of Toeplitz/circulant or introduces the circulant unit shortcut into general Toeplitz.

## Coefficient preparation tradeoff

Both new preparation methods use 3·3ⁿ−4N+1 coefficient XORs. They **do not** retain the old hash layout's coefficient common-subexpression elimination. The old 4×4 counts were 11 for Toeplitz and 5 for circulant; the new direct numeric preparation uses 12 for either. This can cost more coefficient arithmetic while eliminating substantial layout construction, metadata, and allocation overhead. The measured comparison decides the effect; operation counts alone do not.

Both preparation methods produce exactly the same ordered 3ⁿ coefficient bytes as the reference iterative layout. Preparation traversal and application traversal are separate API choices. The benchmark uses recursive preparation with recursive application, and iterative preparation with iterative application.

## API

```cpp
#include "cpp/packed_matrix.hpp"
using namespace dynamic_matrix;

Vector column{2,1,1,3}, x{0xdb,0x13,0x53,0x45};
auto rec = PackedPlan::circulant(column, true, Preparation::recursive);
auto iter = PackedPlan::circulant(column, true, Preparation::iterative);
auto brute_result = brute_force_circulant(column, x);
auto recursive_result = rec.recursive(x);
auto iterative_result = iter.iterative(x);

Vector output(4), work(iter.iterative_workspace_bytes());
iter.iterative_into(x, output, work);  // reuse both for later inputs

auto toeplitz = PackedPlan::toeplitz({2,7,1,3}, {2,9,11,13}, Preparation::iterative);
```

The `_into` APIs require distinct input/output/work vectors and reject insufficient workspace. `PackedPlan` owns its coefficients. Caller buffers must not be shared across simultaneous evaluations without synchronization; separate buffers may reuse a const plan. Existing `Counts` semantics apply and accumulate across calls.

## Costs

Let L=3ⁿ, N=2ⁿ. Counts exclude coefficient preparation and validation.

| Quantity | Brute force | Packed recursive | Packed iterative |
|---|---:|---:|---:|
| General field multiplications | N² | L | L |
| Unit-sum circulant multiplications | N² | L−N | L−N |
| Vector-side XORs | N(N−1) | 3(L−N) | 3(L−N) |
| Retained coefficient bytes | 0 | L | L |
| Workspace excluding output | 0 | 0 for N≤2, else 2N−4 | L |
| Output bytes | N | N | N |
| Calls during application | fixed-depth loops | O(log N) recursive depth | fixed-depth loops |

Byte formulas describe logical/allocated payload in the measured standard-library configuration, excluding allocator metadata, container objects, inputs and stack. The independent profiler checks actual vector capacities and application peaks. Zero per-call allocations in reuse mode do not mean zero memory.

Both preparation paths are O(L) arithmetic/work under the fixed-width field model. Recursive preparation additionally holds 2N−1 input diagonals and scratch of 4N−4−2n bytes for N>1. Iterative preparation holds the L-byte array plus a temporary group of 2N−1 bytes for N>2. No retained size-only hash layout is omitted from the cold timing.

## Validation and reproduction

The C++ suite compares both preparation modes' entire coefficient arrays against the existing layout, and verifies recursive/iterative results, operation counts, changing inputs, zero vectors, input errors, alias rejection, workspace guard bytes, all 72 shared fixtures, and the independent paper formula. Cases cover all three matrix types at powers of two through N=1024. The scalar oracle exhaustively checks all 65,536 byte products.

```sh
c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic cpp/test.cpp -o build/dm_test
./build/dm_test tests/fixtures.txt
c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic cpp/compare_packed.cpp -o build/dm_compare_packed
python tools/run_iterative_comparison.py --profile packed --executable build/dm_compare_packed --iterations 32
c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic tools/packed_memory.cpp -o build/dm_packed_memory
./build/dm_packed_memory > analysis/packed/memory.csv
python tools/report_packed_comparison.py --plots
```

Use `.exe` on Windows. Plot generation requires matplotlib. CMake builds the correctness test and both comparison executables. Local validation used direct Clang commands; remote CI has not been run.

The benchmark uses one fixed matrix and 16 rotating vectors per size. Allocating prepared application, allocation-free prepared application, matrix preparation, and first-use end-to-end costs are separate phases. Old implementations are timed within the same run, so comparisons do not rely on mismatched historical workloads. Raw samples, medians/quartiles, machine/compiler metadata and source hashes are saved in [the new report directory](../analysis/packed/README.md). Historical results remain in `analysis/iterative`.
