#pragma once
// Requested heap payload only; single thread, no over-aligned allocations here.
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
    auto* h = static_cast<Header*>(std::malloc(sizeof(Header)+bytes));
    if (!h) throw std::bad_alloc();
    h->bytes = bytes;
    h->epoch = active ? epoch : 0;
    if (active) { live += bytes; if (live > peak) peak = live; ++allocations; }
    return h+1;
}
void release(void* pointer) noexcept {
    if (!pointer) return;
    auto* h = static_cast<Header*>(pointer)-1;
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
