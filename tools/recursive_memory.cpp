// Measures requested heap bytes, not RSS or allocator metadata. Single-threaded.
#include <cstddef>
#include <cstdlib>
#include <new>

namespace tracked {
struct alignas(std::max_align_t) Header { std::size_t bytes, epoch; };
std::size_t epoch = 0, live = 0, peak = 0, allocations = 0;
bool active = false;
void begin() { ++epoch; live = peak = allocations = 0; active = true; }
void end() { active = false; }
void* allocate(std::size_t bytes) {
    auto* h = static_cast<Header*>(std::malloc(sizeof(Header) + bytes));
    if (!h) throw std::bad_alloc();
    h->bytes = bytes;
    h->epoch = active ? epoch : 0;
    if (active) {
        live += bytes;
        if (live > peak) peak = live;
        ++allocations;
    }
    return h + 1;
}
void release(void* pointer) noexcept {
    if (!pointer) return;
    auto* h = static_cast<Header*>(pointer) - 1;
    if (active && h->epoch == epoch) live -= h->bytes;
    std::free(h);
}
}
void* operator new(std::size_t n) { return tracked::allocate(n); }
void* operator new[](std::size_t n) { return tracked::allocate(n); }
void operator delete(void* p) noexcept { tracked::release(p); }
void operator delete[](void* p) noexcept { tracked::release(p); }
void operator delete(void* p, std::size_t) noexcept { tracked::release(p); }
void operator delete[](void* p, std::size_t) noexcept { tracked::release(p); }

#include "../cpp/dynamic_matrix.hpp"
#include <iostream>
#include <random>
using namespace dynamic_matrix;

int main() {
    std::cout << "n,N,node_count,node_size_bytes,plan_retained_heap_bytes,prepare_peak_heap_bytes,prepare_allocations,apply_peak_extra_heap_bytes,apply_allocations,end_to_end_peak_heap_bytes,baseline_peak_extra_heap_bytes,recursion_depth\n";
    std::mt19937 rng(20261007);
    std::size_t power3 = 1;
    for (int n = 0; n <= 10; ++n, power3 *= 3) {
        std::size_t size = std::size_t(1) << n;
        Vector column(size), x(size);
        for (auto& a : column) a = static_cast<Byte>(rng() & 255);
        for (auto& a : x) a = static_cast<Byte>(rng() & 255);
        tracked::begin();
        auto plan = prepare_circulant(column);
        const auto retained = tracked::live, prepare_peak = tracked::peak, prep_allocs = tracked::allocations;
        tracked::end();
        tracked::begin();
        auto output = plan.matvec(x);
        const auto apply_peak = tracked::peak, apply_allocs = tracked::allocations;
        tracked::end();
        tracked::begin();
        auto one_shot = prepare_circulant(column).matvec(x);
        const auto end_peak = tracked::peak;
        tracked::end();
        tracked::begin();
        auto expected = circulant_matvec(column, x);
        const auto baseline_peak = tracked::peak;
        tracked::end();
        if (output != expected || one_shot != expected) return 1;
        const auto nodes = (3 * power3 - 1) / 2;
        // Validate the profiler against independent allocation-size formulas.
        if (retained != nodes * sizeof(detail::Node) || apply_peak != 3 * size - 2 || baseline_peak != size) return 2;
        std::cout << n << ',' << size << ',' << nodes << ',' << sizeof(detail::Node) << ','
                  << retained << ',' << prepare_peak << ',' << prep_allocs << ',' << apply_peak << ','
                  << apply_allocs << ',' << end_peak << ',' << baseline_peak << ',' << n+1 << '\n';
    }
}
