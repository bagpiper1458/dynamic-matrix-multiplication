#include "dynamic_matrix.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

using namespace dynamic_matrix;
volatile std::uint64_t benchmark_sink = 0;
void consume(const Vector& x) {
    std::uint64_t sum = 0;
    for (Byte a : x) sum += a;
    benchmark_sink = sum;
}
template<class F> double sample(F function, int iterations) {
    const auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < iterations; ++i) function();
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-start).count()/iterations;
}
int integer(const char* value) {
    std::size_t used;
    int result = std::stoi(value, &used);
    if (used != std::string(value).size()) throw std::invalid_argument("invalid integer argument");
    return result;
}
int main(int argc, char** argv) {
    try {
        const int max_n = argc > 1 ? integer(argv[1]) : 7;
        const int repeats = argc > 2 ? integer(argv[2]) : 7;
        const int iterations = argc > 3 ? integer(argv[3]) : 100;
        const std::string kind = argc > 4 ? argv[4] : "circulant";
        if (max_n < 0 || max_n > 10 || repeats < 1 || iterations < 1 || argc > 5 ||
            (kind != "toeplitz" && kind != "circulant" && kind != "unit-circulant"))
            throw std::invalid_argument("usage: dm_benchmark [max_n 0..10] [repeats>0] [iterations>0] [toeplitz|circulant|unit-circulant]");
#if defined(__clang__)
        std::cerr << "compiler=Clang " << __clang_version__ << '\n';
#elif defined(__GNUC__)
        std::cerr << "compiler=GCC " << __VERSION__ << '\n';
#elif defined(_MSC_VER)
        std::cerr << "compiler=MSVC " << _MSC_VER << '\n';
#endif
        std::cerr << "seed=20261007 repeats=" << repeats << " iterations=" << iterations
                  << " scalar_modulus=0x11B units=median_ms counters=disabled\n";
        std::cout << "n,N,kind,baseline_ms,prepare_ms,apply_ms,end_to_end_ms,apply_speedup,end_to_end_speedup\n";
        std::cout << std::fixed << std::setprecision(8);
        std::mt19937 rng(20261007);
        for (int n = 0; n <= max_n; ++n) {
            const std::size_t size = std::size_t(1) << n;
            Vector column(size), row(size), x(size);
            for (auto& a : column) a = static_cast<Byte>(rng() & 255);
            for (auto& a : row) a = static_cast<Byte>(rng() & 255);
            for (auto& a : x) a = static_cast<Byte>(rng() & 255);
            if (kind == "unit-circulant") {
                Byte sum = 1;
                for (std::size_t i = 0; i+1 < size; ++i) sum ^= column[i];
                column.back() = sum;
            }
            row[0] = column[0];
            auto prepare = [&]() { return kind == "toeplitz" ? prepare_toeplitz(column, row) : prepare_circulant(column, kind == "unit-circulant"); };
            auto baseline = [&]() { return kind == "toeplitz" ? toeplitz_matvec(column, row, x) : circulant_matvec(column, x); };
            auto plan = prepare();
            if (plan.matvec(x) != baseline()) throw std::runtime_error("correctness check failed");
            auto run = [&](int which) {
                if (which == 0) consume(baseline());
                if (which == 1) { auto p = prepare(); benchmark_sink = p.preparation_additions(); }
                if (which == 2) consume(plan.matvec(x));
                if (which == 3) consume(prepare().matvec(x));
            };
            for (int which = 0; which < 4; ++which) run(which);
            std::array<std::vector<double>, 4> samples;
            std::array<int, 4> order{0,1,2,3};
            for (int r = 0; r < repeats; ++r) {
                std::shuffle(order.begin(), order.end(), rng);
                for (int which : order) samples[which].push_back(sample([&] { run(which); }, iterations));
            }
            std::array<double, 4> median{};
            for (int i = 0; i < 4; ++i) {
                std::sort(samples[i].begin(), samples[i].end());
                median[i] = (samples[i][(repeats-1)/2] + samples[i][repeats/2])/2;
            }
            std::cout << n << ',' << size << ',' << kind;
            for (auto value : median) std::cout << ',' << value;
            std::cout << ',' << median[0]/median[2] << ',' << median[0]/median[3] << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
