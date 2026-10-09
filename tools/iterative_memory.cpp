#include "tracked_heap.hpp"
#include "../cpp/iterative_matrix.hpp"
#include <iostream>
#include <random>
#include <string>
using namespace dynamic_matrix;
int main() {
    std::cout << "n,N,kind,recursive_plan_bytes,recursive_prepare_peak_bytes,recursive_apply_extra_bytes,layout_retained_bytes,layout_build_peak_bytes,iterative_plan_bytes,iterative_prepare_extra_peak_bytes,iterative_apply_extra_bytes,recursive_total_apply_bytes,iterative_total_apply_without_layout_bytes,iterative_total_apply_with_layout_bytes,recursive_prepare_xors,iterative_prepare_xors,brute_force_peak_bytes,recursive_cold_peak_bytes,iterative_cold_peak_bytes\n";
    std::mt19937 rng(20261009);
    for (const std::string kind : {"toeplitz", "circulant", "unit-circulant"}) {
        std::size_t leaves = 1;
        for (int n = 0; n <= 10; ++n, leaves *= 3) {
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
            tracked::begin();
            const auto expected = circular ? circulant_matvec(column, x) : toeplitz_matvec(column, row, x);
            const auto brute_peak = tracked::peak;
            tracked::end();
            tracked::begin();
            const auto cold_recursive = (circular ? prepare_circulant(column, unit) : prepare_toeplitz(column, row)).matvec(x);
            const auto rec_cold_peak = tracked::peak;
            tracked::end();
            tracked::begin();
            const auto cold_iterative = [&] {
                IterativeLayout cold_layout(size, circular);
                return cold_layout.prepare(column, circular ? empty : row, unit).matvec(x);
            }();
            const auto iter_cold_peak = tracked::peak;
            tracked::end();
            tracked::begin();
            auto recursive = circular ? prepare_circulant(column, unit) : prepare_toeplitz(column, row);
            const auto rec_retained = tracked::live, rec_prep_peak = tracked::peak;
            tracked::end();
            tracked::begin();
            auto r = recursive.matvec(x);
            const auto rec_apply = tracked::peak;
            tracked::end();
            tracked::begin();
            IterativeLayout layout(size, circular);
            const auto layout_retained = tracked::live, layout_peak = tracked::peak;
            tracked::end();
            tracked::begin();
            auto iterative = layout.prepare(column, circular ? empty : row, unit);
            const auto iter_retained = tracked::live, iter_prep_peak = tracked::peak;
            tracked::end();
            tracked::begin();
            auto y = iterative.matvec(x);
            const auto iter_apply = tracked::peak;
            tracked::end();
            if (r != expected || y != expected || cold_recursive != expected || cold_iterative != expected ||
                brute_peak != size || rec_retained != ((3*leaves-1)/2)*sizeof(detail::Node) ||
                rec_apply != 3*size-2 || iter_apply != iterative.application_buffer_bytes() ||
                iter_retained != iterative.retained_payload_bytes() || layout_retained != layout.retained_payload_bytes()) {
                std::cerr << "memory or correctness invariant failed\n";
                return 1;
            }
            std::cout << n << ',' << size << ',' << kind << ',' << rec_retained << ',' << rec_prep_peak << ','
                << rec_apply << ',' << layout_retained << ',' << layout_peak << ',' << iter_retained << ','
                << iter_prep_peak << ',' << iter_apply << ',' << rec_retained+rec_apply << ','
                << iter_retained+iter_apply << ',' << layout_retained+iter_retained+iter_apply << ','
                << recursive.preparation_additions() << ',' << iterative.preparation_additions() << ','
                << brute_peak << ',' << rec_cold_peak << ',' << iter_cold_peak << '\n';
        }
    }
}
