#include "tracked_heap.hpp"
#include "../cpp/packed_matrix.hpp"
#include <iostream>
#include <random>
#include <string>
using namespace dynamic_matrix;
int main() {
    std::cout << "exponent,N,kind,brute_force_peak_bytes,recursive_plan_bytes,iterative_plan_bytes,recursive_prepare_peak_bytes,iterative_prepare_peak_bytes,recursive_apply_extra_bytes,iterative_apply_extra_bytes,recursive_total_apply_bytes,iterative_total_apply_bytes,recursive_cold_peak_bytes,iterative_cold_peak_bytes,brute_reuse_allocations,recursive_reuse_allocations,iterative_reuse_allocations\n";
    std::mt19937 rng(20261009);
    for (const std::string kind : {"toeplitz", "circulant", "unit-circulant"}) {
        for (int exponent = 0; exponent <= 10; ++exponent) {
            const auto n = std::size_t(1) << exponent;
            const bool circular = kind != "toeplitz", unit = kind == "unit-circulant";
            Vector column(n), row(n), x(n);
            for (auto& a : column) a = static_cast<Byte>(rng() & 255);
            for (auto& a : row) a = static_cast<Byte>(rng() & 255);
            for (auto& a : x) a = static_cast<Byte>(rng() & 255);
            if (unit) { Byte sum = 1; for (std::size_t i = 0; i+1 < n; ++i) sum ^= column[i]; column.back() = sum; }
            row[0] = column[0];
            auto prepare = [&](Preparation method) {
                return circular ? PackedPlan::circulant(column,unit,method) : PackedPlan::toeplitz(column,row,method);
            };
            tracked::begin();
            auto truth = circular ? brute_force_circulant(column,x) : brute_force_toeplitz(column,row,x);
            auto bf_peak = tracked::peak;
            tracked::end();
            tracked::begin(); auto recursive = prepare(Preparation::recursive);
            auto rp = tracked::peak, rr = tracked::live; tracked::end();
            tracked::begin(); auto iterative = prepare(Preparation::iterative);
            auto ip = tracked::peak, ir = tracked::live; tracked::end();
            tracked::begin(); auto r = recursive.recursive(x);
            auto ra = tracked::peak; tracked::end();
            tracked::begin(); auto y = iterative.iterative(x);
            auto ia = tracked::peak; tracked::end();
            tracked::begin(); auto cold_r = prepare(Preparation::recursive).recursive(x);
            auto rcold = tracked::peak; tracked::end();
            tracked::begin(); auto cold_i = prepare(Preparation::iterative).iterative(x);
            auto icold = tracked::peak; tracked::end();
            Vector out(n), rw(recursive.recursive_workspace_bytes()), iw(iterative.iterative_workspace_bytes());
            tracked::begin();
            if (circular) brute_force_circulant_into(column,x,out); else brute_force_toeplitz_into(column,row,x,out);
            auto ba_count = tracked::allocations; tracked::end();
            if (out != truth) return 1;
            tracked::begin(); recursive.recursive_into(x,out,rw);
            auto ra_count = tracked::allocations; tracked::end();
            if (out != truth) return 2;
            tracked::begin(); iterative.iterative_into(x,out,iw);
            auto ia_count = tracked::allocations; tracked::end();
            if (out != truth || r != truth || y != truth || cold_r != truth || cold_i != truth ||
                rr != recursive.retained_payload_bytes() || ir != iterative.retained_payload_bytes() ||
                rr != recursive.coefficients().size() || ir != rr || bf_peak != n ||
                ra != n+rw.size() || ia != n+iw.size() || ba_count || ra_count || ia_count) return 3;
            std::cout << exponent << ',' << n << ',' << kind << ',' << bf_peak << ',' << rr << ',' << ir << ','
                      << rp << ',' << ip << ',' << ra << ',' << ia << ',' << rr+ra << ',' << ir+ia << ','
                      << rcold << ',' << icold << ',' << ba_count << ',' << ra_count << ',' << ia_count << '\n';
        }
    }
}
