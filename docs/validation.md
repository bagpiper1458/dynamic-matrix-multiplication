# Local validation record

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
