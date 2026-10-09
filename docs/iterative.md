# Iterative PARAMS reconstruction

The author reviewed the [4×4 formula/schedule audit](../analysis/params4/index.html) before this C++ implementation was generalized. The author later supplied recovered formula-generation code and clarified that PARAMS reuse had only existed as pseudocode. This is a new implementation of those ideas and does not claim a globally optimal schedule.

See [paper-to-code correspondence (Traditional Chinese)](paper-correspondence.zh-TW.md) for the manuscript's SM_A/SM_B generation, the F/G/H identities, and independent symbolic checks. This implementation compacts the same formulas into a ternary schedule; it is not a line-by-line translation of Algorithms 1–4.

The subsequent [PARAMS derivation](params-derivation.zh-TW.md) defines aligned block values and exact structural keys. Its child-key relation (3t, 3t+1, 3t+2) explains how the existing contraction loop evaluates the derived PARAMS blocks. This establishes a concrete correspondence without asserting that the manuscript uniquely specified the storage layout.

Both Python (`dynamic_matrix/iterative.py`) and C++ (`cpp/iterative_matrix.hpp`) have three separate stages, with no recursive calls in their iterative paths:

1. **Layout, once per size and matrix structure.** Generate a coefficient DAG using canonical operand-ID pairs, plus ternary product order and the structural unit-leaf indices. The temporary hash table is discarded after construction. No input coefficient values participate in expression sharing.
2. **Prepare, once per matrix.** Run the coefficient XOR DAG and gather the 3ⁿ product coefficients into a byte array. Validate the unit-row-sum condition when requested. The resulting plan owns its data and can outlive its layout.
3. **Apply, once per vector.** Expand binary inputs backwards into shared subcube sums, multiply corresponding coefficients, then contract ternary products forwards into outputs. Each in-place step reads its operands before overwriting them. Only the structurally justified unit leaves skip multiplication.

The 4×4 plan has nine products and 15 vector-side XORs. General Toeplitz coefficient preparation uses 11 XORs; circulant preparation uses five. A unit-sum circulant performs five field multiplications and four identity copies. Validation is separate from these arithmetic counts.

## APIs

```python
from dynamic_matrix import make_iterative_layout

layout = make_iterative_layout(4, circulant=True)
plan = layout.prepare([2, 1, 1, 3], unit_row_sum=True)
assert plan.matvec([0xDB, 0x13, 0x53, 0x45]) == (0x8E, 0x4D, 0xA1, 0xBC)
# Same layout, different matrix:
assert layout.prepare([1, 0, 0, 0]).matvec([1, 2, 3, 4]) == (1, 2, 3, 4)
```

```cpp
#include "cpp/iterative_matrix.hpp"
using namespace dynamic_matrix;

IterativeLayout circular(4, true);
auto plan = circular.prepare({2, 1, 1, 3}, {}, true);
auto result = plan.matvec({0xDB, 0x13, 0x53, 0x45});

IterativeLayout toeplitz(4);
auto other = toeplitz.prepare({2, 7, 1, 3}, {2, 9, 11, 13});
```

C++ parameters are `(first_column, first_row = {}, unit_row_sum = false)`. A circulant layout requires an empty row; a Toeplitz layout requires a matching row and does not accept the unit shortcut. Both APIs check lengths. C++ uses byte types; validate wider numbers before narrowing.

## Complexity and memory boundaries

Let N=2ⁿ, L=3ⁿ, Q be the coefficient DAG's XOR count, and U=N for a unit-sum circulant or zero otherwise.

| Item | Recursive C++ | Iterative C++ |
|---|---|---|
| Field multiplications per vector | L−U | L−U |
| Vector-side XORs | 3(L−N) | 3(L−N) |
| Matrix preparation XORs | 3L−4N+1 | Q |
| Prepared matrix storage | (3L−1)/2 heap nodes | L coefficient bytes and U `size_t` indices |
| Application heap buffers, including output | 3N−2 bytes | L+N bytes |
| Algorithmic call depth | n+1 | constant |
| Separate reusable layout | none | Q operand-ID pairs, L leaf IDs, N unit-leaf IDs |

Counts exclude input checks and integer indexing; storage formulas exclude allocator overhead, stack frames, and container objects. The profiler includes actual retained vector **capacities**, which can exceed logical lengths. It reports layout build peak separately, including temporary hash-table and stage buffers.

Application time is O(L) with a fixed eight-step field kernel in both versions. Iteration alone does not reduce the multiplication or addition order. Layout construction traverses O(L) entries; hash interning has expected O(1) lookup, not a worst-case constant-time guarantee. Matrix preparation and gathering cost O(N+Q+L). Layout and plan storage are O(L), as is this iterative work buffer. The recursive evaluator has O(N) extra buffers and O(log N) call depth, so it has a smaller **per-call extra workspace** at sufficiently large N. The iterative plan's compact coefficient array can nevertheless substantially reduce **total retained-plus-working memory** compared with the pointer tree.

Sharing canonical pairs does not detect every algebraically equivalent XOR expression. This version retains duplicate leaf coefficient bytes for contiguous evaluation, and does not claim minimal memory or the fastest possible implementation. The prepared plan is already independent of the layout: discard the layout when no more matrices need preparation, or share it across multiple plans of the same structure and size.

## Reproduce the comparison

```sh
c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic cpp/compare.cpp -o build/dm_compare
python tools/run_iterative_comparison.py --executable build/dm_compare
c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic tools/iterative_memory.cpp -o build/dm_iterative_memory
./build/dm_iterative_memory > analysis/iterative/memory.csv
python tools/report_iterative_comparison.py
# Optional PNG/SVG figures and standalone HTML (requires matplotlib):
python tools/report_iterative_comparison.py --plots
```

On Windows use `.exe` paths. The timing runner saves all samples and source SHA-256 fingerprints. Timing and heap profiling are separate executables. See the [recorded report](../analysis/iterative/README.md) for measurement boundaries and this machine's results.
