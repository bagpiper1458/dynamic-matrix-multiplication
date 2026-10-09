# Local validation record

## Publication checks — 2026-10-09 (Asia/Taipei)

- Python 3.12.14: all 17 tests passed again before upload; all 72 checked-in fixtures match regeneration by the independent oracle.
- Published benchmark metadata uses repository-relative executable paths; timing samples, run arguments, compiler details and source hashes are unchanged.
- The GitHub Actions badge reports remote status separately from the local validation records below.

## Compact C++ implementations — 2026-10-09 (Asia/Taipei)

- C++17, Clang 18.1.6 via Zig 0.13.0, `-O3 -Wall -Wextra -Wpedantic`: the updated suite passed after the final inner-loop optimization.
- Both compact preparation methods produce the same complete ordered coefficient arrays as the reference iterative layout. Both compact evaluators and improved brute force match independent/reference results for Toeplitz, circulant, and unit-sum circulant at powers of two through N=1024.
- Checks cover all 65,536 scalar byte products, 72 shared fixtures, the independent paper formula through N=16, field operation counts, changed and zero inputs, alias rejection, invalid sizes, insufficient workspace, and workspace guard bytes.
- The independent heap profiler was rebuilt from the final header and passed all three matrix types through N=1024. Every reuse API made zero allocations. At N=1024, each plan retains 59,049 coefficient bytes; prepared recursive evaluation additionally uses 3,068 bytes including output, and iterative evaluation uses 60,073. These are requested heap payloads, not RSS; input, stack, container objects and allocator metadata are excluded.
- Final timing used seven repeats of 32 calls per phase, one fixed matrix and 16 rotating vectors, at sizes 1–1024. The old and new implementations were measured in the same run. Heap profiling and compilation were not run alongside timing. [Raw samples, metadata/source hashes, summary quartiles and charts](../analysis/packed/README.md) preserve the result.
- Source hashes match the final measured C++ sources; the archive contains 3,465 unique timing samples and 33 memory cases. Both static chart images were visually inspected. The generated HTML embeds these charts; no automated browser verification is claimed.
- Python algorithms were unchanged in this optimization; the preceding Python validation below remains applicable. Remote CI has not been run.

## Iterative extension — 2026-10-09 (Asia/Taipei)

- Python 3.12.14: all 17 tests passed, including literal general-formula comparisons, independent dense arithmetic, 4×4 symbolic cancellation, and register-lifetime replay. The additional manuscript Algorithm 1/2 tests generate sets independently, validate symbolic coefficient correspondence and bilinear cancellation through N=32, and numerically check both paths through N=16.
- C++17, Clang 18.1.6 via Zig 0.13.0, `-O3 -Wall -Wextra -Wpedantic`: recursive and iterative results and operation counts matched for all three matrix types at N=1,2,...,1024 (powers of two).
- Both C++ paths passed all 72 independently generated shared fixtures and all 65,536 scalar products. Tests cover layout reuse with a new matrix, plan use after layout destruction, zero/identity matrices, malformed inputs, and the reviewed 4×4 preparation counts (11/5).
- C++ also independently generates Algorithm 1/2's sets and checks their numeric formula through N=16; the approved 4×4 coefficient ordering is checked explicitly.
- Separate heap profiler checked all three matrix types through N=1024. It asserts measured retained storage against vector-capacity/node formulas and application peaks against 3N−2 and 3ⁿ+N.
- The revised three-way report includes brute-force heap allocation (N bytes, including output) and measured complete cold-start peaks for recursive and iterative evaluation. Every chart compares all three algorithms under the same stated boundary. Timing samples are unchanged; additional memory phases were measured separately.
- Fresh timing runs used seven repeats of 20 calls per phase, sizes 1–1024, all three matrix types. Raw samples, compiler/platform metadata, source hashes, summary quartiles, and PNG/SVG charts are in [analysis/iterative](../analysis/iterative/README.md). Static charts were visually inspected; no automated browser preview is claimed.

These are local results. Remote CI has not been run. The initial recursive measurements below are preserved as a historical record.

## Initial recursive reconstruction

Verified on 2026-10-07 (Asia/Taipei), Windows 11 x86-64.

| Check | Local result |
|---|---|
| Python 3.12.14 `unittest` | 9 tests passed; tests contain exhaustive and parameterized checks |
| Python finite-field oracle | All 65,536 byte products matched independent polynomial long division |
| Python matrix tests | Toeplitz/general circulant/unit-sum circulant through N=128 matched a dense polynomial-arithmetic oracle |
| Python/C++ shared fixtures | 72 deterministic cases, sizes 1 through 128 |
| C++17 compilation | Clang 18.1.6 via Zig 0.13.0; `-O3 -Wall -Wextra -Wpedantic` |
| C++ test executable | Passed 65,536 scalar products, matrix/count/invalid-input tests, and 72 shared fixtures |
| Operation counts | General 3ⁿ M, unit 3ⁿ−2ⁿ M, both 3(3ⁿ−2ⁿ) A; tested n=0 through 7 |
| Benchmarks | Both languages and all three matrix variants; correctness checked before timing |

The first compiler invocation emitted warnings from Zig's bundled unwind library, not from project source files. Compilation and test execution succeeded.

Reproducible commands:

```sh
python -m unittest discover -s tests -v
python tools/generate_fixtures.py --check
c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic cpp/test.cpp -o dm_test
./dm_test tests/fixtures.txt
```

The CMake and GitHub Actions configurations are supplied for other environments; the local C++ verification used direct compiler commands. A local pass does not establish that a remote Actions run or all other compilers have passed. Check the repository's Actions page after publishing.
