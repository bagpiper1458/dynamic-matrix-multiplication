#include "iterative_matrix.hpp"
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#endif

using namespace dynamic_matrix;
volatile std::uint64_t comparison_sink = 0;
const void* volatile comparison_pointer = nullptr;
void escape(const void* pointer) {
#if defined(__clang__) || defined(__GNUC__)
    asm volatile("" : : "g"(pointer) : "memory");
#else
    comparison_pointer = pointer;
    _ReadWriteBarrier();
#endif
}
void consume(const Vector& x) {
    std::uint64_t sum = 0;
    for (Byte a : x) sum += a;
    comparison_sink = sum;
}
int integer(const char* value) {
    std::size_t used;
    const int result = std::stoi(value, &used);
    if (used != std::string(value).size()) throw std::invalid_argument("invalid integer argument");
    return result;
}
int main(int argc, char** argv) {
    try {
        const int max_n = argc > 1 ? integer(argv[1]) : 10;
        const int repeats = argc > 2 ? integer(argv[2]) : 7;
        const int iterations = argc > 3 ? integer(argv[3]) : 20;
        const std::string kind = argc > 4 ? argv[4] : "circulant";
        if (argc > 5 || max_n < 0 || max_n > 10 || repeats < 1 || iterations < 1 ||
            (kind != "toeplitz" && kind != "circulant" && kind != "unit-circulant"))
            throw std::invalid_argument("usage: dm_compare [max_n 0..10] [repeats>0] [iterations>0] [toeplitz|circulant|unit-circulant]");
#if defined(__clang__)
        std::cerr << "compiler=Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
        std::cerr << "compiler=GCC " << __VERSION__ << '\n';
#elif defined(_MSC_VER)
        std::cerr << "compiler=MSVC " << _MSC_VER << '\n';
#endif
        std::cerr << "seed=20261009 repeats=" << repeats << " iterations=" << iterations
                  << " scalar_modulus=0x11B units=ms_per_call counters=disabled sizeof_size_t=" << sizeof(std::size_t)
                  << " sizeof_recursive_node=" << sizeof(detail::Node) << '\n';
        const std::array<const char*, 9> names{{"baseline", "recursive_prepare", "recursive_apply", "recursive_end",
            "iterative_layout", "iterative_prepare", "iterative_apply", "iterative_hot_end", "iterative_cold_end"}};
        std::cout << "n,N,kind,repeat,phase,ms_per_call\n" << std::fixed << std::setprecision(9);
        std::mt19937 rng(20261009);
        for (int n = 0; n <= max_n; ++n) {
            const std::size_t size = std::size_t(1) << n;
            Vector column(size), row(size), x(size);
            for (auto& a : column) a = static_cast<Byte>(rng() & 255);
            for (auto& a : row) a = static_cast<Byte>(rng() & 255);
            for (auto& a : x) a = static_cast<Byte>(rng() & 255);
            const bool circular = kind != "toeplitz", unit = kind == "unit-circulant";
            if (unit) {
                Byte sum = 1;
                for (std::size_t i = 0; i+1 < size; ++i) sum ^= column[i];
                column.back() = sum;
            }
            row[0] = column[0];
            const Vector empty;
            auto recursive_prepare = [&] { return circular ? prepare_circulant(column, unit) : prepare_toeplitz(column, row); };
            auto baseline = [&] { return circular ? circulant_matvec(column, x) : toeplitz_matvec(column, row, x); };
            const IterativeLayout layout(size, circular);
            auto iterative_prepare = [&] { return layout.prepare(column, circular ? empty : row, unit); };
            const auto recursive = recursive_prepare();
            const auto iterative = iterative_prepare();
            const auto expected = baseline();
            if (recursive.matvec(x) != expected || iterative.matvec(x) != expected)
                throw std::runtime_error("pre-timing correctness check failed");
            auto run = [&](int which) {
                switch (which) {
                case 0: consume(baseline()); break;
                case 1: { auto p = recursive_prepare(); escape(&p); break; }
                case 2: consume(recursive.matvec(x)); break;
                case 3: consume(recursive_prepare().matvec(x)); break;
                case 4: { IterativeLayout l(size, circular); escape(&l); break; }
                case 5: { auto p = iterative_prepare(); escape(&p); break; }
                case 6: consume(iterative.matvec(x)); break;
                case 7: consume(iterative_prepare().matvec(x)); break;
                case 8: { IterativeLayout l(size, circular); consume(l.prepare(column, circular ? empty : row, unit).matvec(x)); break; }
                }
            };
            std::array<int, 9> order{{0,1,2,3,4,5,6,7,8}};
            for (auto which : order) run(which);
            for (int r = 0; r < repeats; ++r) {
                std::shuffle(order.begin(), order.end(), rng);
                for (auto which : order) {
                    const auto start = std::chrono::steady_clock::now();
                    for (int i = 0; i < iterations; ++i) run(which);
                    const auto ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-start).count()/iterations;
                    std::cout << n << ',' << size << ',' << kind << ',' << r << ',' << names[which] << ',' << ms << '\n';
                }
            }
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
