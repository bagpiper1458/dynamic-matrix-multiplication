#include "packed_matrix.hpp"
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
    int result = std::stoi(value, &used);
    if (used != std::string(value).size()) throw std::invalid_argument("invalid integer");
    return result;
}
int main(int argc, char** argv) {
    try {
        const int maximum = argc > 1 ? integer(argv[1]) : 10;
        const int repeats = argc > 2 ? integer(argv[2]) : 7;
        const int iterations = argc > 3 ? integer(argv[3]) : 32;
        const std::string kind = argc > 4 ? argv[4] : "circulant";
        if (argc > 5 || maximum < 0 || maximum > 10 || repeats < 1 || iterations < 1 ||
            (kind != "toeplitz" && kind != "circulant" && kind != "unit-circulant"))
            throw std::invalid_argument("usage: dm_compare_packed [max_n 0..10] [repeats>0] [iterations>0] [kind]");
#if defined(__clang__)
        std::cerr << "compiler=Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
        std::cerr << "compiler=GCC " << __VERSION__ << '\n';
#elif defined(_MSC_VER)
        std::cerr << "compiler=MSVC " << _MSC_VER << '\n';
#endif
        std::cerr << "seed=20261009 repeats=" << repeats << " iterations=" << iterations
                  << " vectors_per_size=16 units=ms_per_call counters=disabled field_kernel=gf_mul_0x11B\n";
        const std::array<const char*, 15> names{{"baseline", "recursive_prepare", "iterative_prepare",
            "recursive_apply", "iterative_apply", "recursive_end", "iterative_end",
            "baseline_reuse", "recursive_reuse", "iterative_reuse", "reference_baseline",
            "reference_recursive_apply", "reference_iterative_apply", "reference_recursive_end", "reference_iterative_cold"}};
        std::mt19937 rng(20261009);
        std::cout << "exponent,N,kind,repeat,phase,ms_per_call\n" << std::fixed << std::setprecision(9);
        for (int exponent = 0; exponent <= maximum; ++exponent) {
            const auto n = std::size_t(1) << exponent;
            const bool circular = kind != "toeplitz", unit = kind == "unit-circulant";
            Vector column(n), row(n), empty;
            for (auto& a : column) a = static_cast<Byte>(rng() & 255);
            for (auto& a : row) a = static_cast<Byte>(rng() & 255);
            if (unit) { Byte sum = 1; for (std::size_t i = 0; i+1 < n; ++i) sum ^= column[i]; column.back() = sum; }
            row[0] = column[0];
            std::array<Vector, 16> inputs;
            for (auto& x : inputs) { x.resize(n); for (auto& a : x) a = static_cast<Byte>(rng() & 255); }
            auto prepare = [&](Preparation method) {
                return circular ? PackedPlan::circulant(column, unit, method) : PackedPlan::toeplitz(column, row, method);
            };
            auto old_prepare = [&] { return circular ? prepare_circulant(column, unit) : prepare_toeplitz(column, row); };
            const auto recursive = prepare(Preparation::recursive), iterative = prepare(Preparation::iterative);
            const auto old_recursive = old_prepare();
            const IterativeLayout old_layout(n, circular);
            const auto old_iterative = old_layout.prepare(column, circular ? empty : row, unit);
            Vector out(n), rec_work(recursive.recursive_workspace_bytes()), iter_work(iterative.iterative_workspace_bytes());
            auto brute = [&](const Vector& x) { return circular ? brute_force_circulant(column,x) : brute_force_toeplitz(column,row,x); };
            auto brute_into = [&](const Vector& x) {
                if (circular) brute_force_circulant_into(column,x,out); else brute_force_toeplitz_into(column,row,x,out);
            };
            auto old_brute = [&](const Vector& x) { return circular ? circulant_matvec(column,x) : toeplitz_matvec(column,row,x); };
            if (recursive.coefficients() != iterative.coefficients() || iterative.coefficients() != old_iterative.coefficients())
                throw std::runtime_error("coefficient ordering mismatch before timing");
            for (const auto& x : inputs) {
                const auto expected = old_brute(x);
                if (brute(x) != expected || recursive.recursive(x) != expected || iterative.iterative(x) != expected ||
                    old_recursive.matvec(x) != expected || old_iterative.matvec(x) != expected)
                    throw std::runtime_error("result mismatch before timing");
                brute_into(x); if (out != expected) throw std::runtime_error("brute reuse mismatch");
                recursive.recursive_into(x,out,rec_work); if (out != expected) throw std::runtime_error("recursive reuse mismatch");
                iterative.iterative_into(x,out,iter_work); if (out != expected) throw std::runtime_error("iterative reuse mismatch");
            }
            auto run = [&](int phase, std::size_t index) {
                const auto& x = inputs[index%inputs.size()];
                switch (phase) {
                case 0: consume(brute(x)); break;
                case 1: { auto p = prepare(Preparation::recursive); escape(&p); break; }
                case 2: { auto p = prepare(Preparation::iterative); escape(&p); break; }
                case 3: consume(recursive.recursive(x)); break;
                case 4: consume(iterative.iterative(x)); break;
                case 5: consume(prepare(Preparation::recursive).recursive(x)); break;
                case 6: consume(prepare(Preparation::iterative).iterative(x)); break;
                case 7: brute_into(x); consume(out); break;
                case 8: recursive.recursive_into(x,out,rec_work); consume(out); break;
                case 9: iterative.iterative_into(x,out,iter_work); consume(out); break;
                case 10: consume(old_brute(x)); break;
                case 11: consume(old_recursive.matvec(x)); break;
                case 12: consume(old_iterative.matvec(x)); break;
                case 13: consume(old_prepare().matvec(x)); break;
                case 14: { IterativeLayout l(n,circular); consume(l.prepare(column,circular ? empty : row,unit).matvec(x)); break; }
                }
            };
            std::array<int, 15> order{};
            for (int i = 0; i < 15; ++i) { order[i] = i; run(i,0); }
            for (int r = 0; r < repeats; ++r) {
                std::shuffle(order.begin(), order.end(), rng);
                for (auto phase : order) {
                    const auto start = std::chrono::steady_clock::now();
                    for (int i = 0; i < iterations; ++i) run(phase, static_cast<std::size_t>(r)*iterations+i);
                    const auto ms = std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count()/iterations;
                    std::cout << exponent << ',' << n << ',' << kind << ',' << r << ',' << names[phase] << ',' << ms << '\n';
                }
            }
        }
    } catch (const std::exception& error) { std::cerr << error.what() << '\n'; return 1; }
}
