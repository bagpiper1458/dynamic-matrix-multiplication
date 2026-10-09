#pragma once

// Compact implementations sharing the scalar field kernel and leaf ordering.
// The older pointer-tree and hash-layout implementations remain as references.
#include "dynamic_matrix.hpp"
#include <algorithm>
#include <limits>

namespace dynamic_matrix {
enum class Preparation { recursive, iterative };

namespace packed_detail {
inline std::size_t leaf_count(std::size_t n) {
    detail::check_size(n);
    std::size_t leaves = 1;
    for (auto width = n; width > 1; width /= 2) {
        if (leaves > std::numeric_limits<std::size_t>::max()/64/3)
            throw std::length_error("matrix size overflows coefficient storage");
        leaves *= 3;
    }
    return leaves;
}
inline void separate(const Vector& x, const Vector& out, const Vector& work) {
    if (&x == &out || &x == &work || &out == &work)
        throw std::invalid_argument("input, output, and workspace must be distinct vectors");
}
inline void compile_recursive(const Byte* d, std::size_t n, Byte* coefficients,
                              std::size_t leaves, Byte* scratch) {
    if (n == 1) { coefficients[0] = d[0]; return; }
    const auto width = n-1, half = n/2, children = leaves/3;
    Byte* g = scratch;
    Byte* h = scratch+width;
    Byte* next = scratch+2*width;
    for (std::size_t i = 0; i < width; ++i) {
        g[i] = d[half+i] ^ d[i];
        h[i] = d[half+i] ^ d[n+i];
    }
    compile_recursive(d+half, half, coefficients, children, next);
    compile_recursive(g, half, coefficients+children, children, next);
    compile_recursive(h, half, coefficients+2*children, children, next);
}
inline void compile_iterative(Vector& coefficients, std::size_t n) {
    if (n == 1) return;
    Vector source(n > 2 ? 2*n-1 : 0);
    Byte* data = coefficients.data();
    Byte* temporary = source.data();
    for (std::size_t width = n, groups = 1; width > 1; width /= 2, groups *= 3) {
        const auto child = width-1, half = width/2;
        // Destination regions grow: work backwards, copying this group's
        // diagonal values before overwriting any of them.
        for (std::size_t group = groups; group-- > 0;) {
            const auto from = group*(2*width-1), to = group*3*child;
            if (width == 2) {
                const Byte a = data[from+1], b = data[from], c = data[from+2];
                data[to] = a;
                data[to+1] = a ^ b;
                data[to+2] = a ^ c;
            } else {
                std::copy_n(data+from, 2*width-1, temporary);
                for (std::size_t i = 0; i < child; ++i) {
                    const auto a = temporary[half+i];
                    data[to+i] = a;
                    data[to+child+i] = a ^ temporary[i];
                    data[to+2*child+i] = a ^ temporary[width+i];
                }
            }
        }
    }
}

template<bool Unit>
inline void evaluate_recursive(const Byte* coefficients, std::size_t leaves,
                               const Byte* x, Byte* out, Byte* scratch, std::size_t n) {
    if (n == 1) {
        if constexpr (Unit) out[0] = x[0];
        else out[0] = gf_mul(coefficients[0], x[0]);
        return;
    }
    // Fuse the final 2x2 block to avoid three separate scalar calls.
    if (n == 2) {
        const auto f = gf_mul(coefficients[0], x[0] ^ x[1]);
        if constexpr (Unit) { out[0] = f ^ x[1]; out[1] = f ^ x[0]; }
        else { out[0] = f ^ gf_mul(coefficients[1], x[1]); out[1] = f ^ gf_mul(coefficients[2], x[0]); }
        return;
    }
    const auto half = n/2, child = leaves/3;
    Byte* mixed = scratch;
    Byte* f = scratch+half;
    Byte* next = scratch+n;
    for (std::size_t i = 0; i < half; ++i) mixed[i] = x[i] ^ x[half+i];
    evaluate_recursive<false>(coefficients, child, mixed, f, next, half);
    evaluate_recursive<Unit>(coefficients+child, child, x+half, out, next, half);
    evaluate_recursive<Unit>(coefficients+2*child, child, x, out+half, next, half);
    for (std::size_t i = 0; i < half; ++i) { out[i] ^= f[i]; out[half+i] ^= f[i]; }
}

inline void multiply_range(const Vector& coefficients, Vector& work, std::size_t begin, std::size_t end) {
    // No per-element identity branch; the compiler may auto-vectorize this loop.
    const Byte* source = coefficients.data();
    Byte* target = work.data();
    for (std::size_t i = begin; i < end; ++i) target[i] = gf_mul(source[i], target[i]);
}
inline void evaluate_iterative(const Vector& coefficients, const Vector& x, Vector& out,
                               Vector& work, bool unit) {
    std::copy(x.begin(), x.end(), work.begin());
    Byte* buffer = work.data();
    for (std::size_t stride = 1, groups = x.size()/2; groups; stride *= 3, groups /= 2) {
        for (std::size_t group = groups; group-- > 0;) {
            const auto from = group*2*stride, to = group*3*stride;
            for (std::size_t i = stride; i-- > 0;) {
                const Byte a = buffer[from+i], b = buffer[from+stride+i];
                buffer[to+i] = a ^ b;
                buffer[to+stride+i] = b;
                buffer[to+2*stride+i] = a;
            }
        }
    }
    if (!unit) multiply_range(coefficients, work, 0, coefficients.size());
    else {
        // Unit leaves have ternary digits exclusively 1/2. Generate their
        // sorted positions with a binary counter, without retaining an ID list.
        std::size_t begin = 0, identity = (coefficients.size()-1)/2;
        for (std::size_t ordinal = 0; ordinal < x.size(); ++ordinal) {
            multiply_range(coefficients, work, begin, identity);
            begin = identity+1;
            if (ordinal+1 < x.size()) {
                std::size_t power3 = 1;
                for (auto v = ordinal+1; !(v & 1); v >>= 1) power3 *= 3;
                identity += (power3+1)/2;
            }
        }
        multiply_range(coefficients, work, begin, coefficients.size());
    }
    for (std::size_t stride = 1, groups = coefficients.size()/3; groups; stride *= 2, groups /= 3) {
        for (std::size_t group = 0; group < groups; ++group) {
            const auto from = group*3*stride, to = group*2*stride;
            for (std::size_t i = 0; i < stride; ++i) {
                const Byte a = buffer[from+i], b = buffer[from+stride+i], c = buffer[from+2*stride+i];
                buffer[to+i] = a ^ b;
                buffer[to+stride+i] = a ^ c;
            }
        }
    }
    std::copy_n(buffer, x.size(), out.begin());
}
} // namespace packed_detail

class PackedPlan {
    std::size_t size_, leaves_;
    bool unit_;
    Vector coefficients_;
    PackedPlan(const Vector& column, const Vector& row, bool circular, bool unit, Preparation method)
        : size_(column.size()), leaves_(packed_detail::leaf_count(size_)), unit_(unit), coefficients_(leaves_) {
        if (!circular) detail::check_toeplitz(column, row);
        if (unit) {
            Byte sum = 0;
            for (Byte value : column) sum ^= value;
            if (sum != 1) throw std::invalid_argument("column XOR must equal one");
        }
        if (size_ == 1) { coefficients_[0] = column[0]; return; }
        auto diagonals = [&](Byte* d) {
            for (std::size_t i = 0; i < size_-1; ++i) d[i] = circular ? column[i+1] : row[size_-1-i];
            std::copy(column.begin(), column.end(), d+size_-1);
        };
        if (method == Preparation::recursive) {
            Vector d(2*size_-1);
            diagonals(d.data());
            std::size_t scratch_bytes = 0;
            for (auto width = size_; width > 1; width /= 2) scratch_bytes += 2*(width-1);
            Vector scratch(scratch_bytes);
            packed_detail::compile_recursive(d.data(), size_, coefficients_.data(), leaves_, scratch.data());
        } else {
            diagonals(coefficients_.data());
            packed_detail::compile_iterative(coefficients_, size_);
        }
    }
    void check(const Vector& x, const Vector& out, const Vector& work, std::size_t required) const {
        if (x.size() != size_ || out.size() != size_ || work.size() < required)
            throw std::invalid_argument("input/output size or workspace size mismatch");
        packed_detail::separate(x, out, work);
    }
    void count(Counts* counts) const {
        if (counts) {
            counts->multiplications += leaves_-(unit_ ? size_ : 0);
            counts->additions += 3*(leaves_-size_);
            counts->identity_copies += unit_ ? size_ : 0;
        }
    }
public:
    static PackedPlan toeplitz(const Vector& column, const Vector& row, Preparation method = Preparation::iterative) {
        detail::check_toeplitz(column, row);
        return PackedPlan(column, row, false, false, method);
    }
    static PackedPlan circulant(const Vector& column, bool unit = false, Preparation method = Preparation::iterative) {
        detail::check_size(column.size());
        return PackedPlan(column, {}, true, unit, method);
    }
    std::size_t size() const { return size_; }
    std::size_t preparation_additions() const { return 3*leaves_-4*size_+1; }
    const Vector& coefficients() const { return coefficients_; }
    std::size_t retained_payload_bytes() const { return coefficients_.capacity(); }
    // The fused 2x2 leaf needs no scratch. Sum N+N/2+...+4 = 2N-4.
    std::size_t recursive_workspace_bytes() const { return size_ <= 2 ? 0 : 2*size_-4; }
    std::size_t iterative_workspace_bytes() const { return leaves_; }
    void recursive_into(const Vector& x, Vector& out, Vector& work, Counts* counts = nullptr) const {
        check(x, out, work, recursive_workspace_bytes());
        if (unit_) packed_detail::evaluate_recursive<true>(coefficients_.data(), leaves_, x.data(), out.data(), work.data(), size_);
        else packed_detail::evaluate_recursive<false>(coefficients_.data(), leaves_, x.data(), out.data(), work.data(), size_);
        count(counts);
    }
    void iterative_into(const Vector& x, Vector& out, Vector& work, Counts* counts = nullptr) const {
        check(x, out, work, iterative_workspace_bytes());
        packed_detail::evaluate_iterative(coefficients_, x, out, work, unit_);
        count(counts);
    }
    Vector recursive(const Vector& x, Counts* counts = nullptr) const {
        Vector out(size_), work(recursive_workspace_bytes());
        recursive_into(x, out, work, counts);
        return out;
    }
    Vector iterative(const Vector& x, Counts* counts = nullptr) const {
        Vector out(size_), work(iterative_workspace_bytes());
        iterative_into(x, out, work, counts);
        return out;
    }
};

inline void brute_force_toeplitz_into(const Vector& column, const Vector& row, const Vector& x,
                                     Vector& out, Counts* counts = nullptr) {
    detail::check_toeplitz(column, row);
    const auto n = column.size();
    if (x.size() != n || out.size() != n || &out == &x || &out == &column || &out == &row)
        throw std::invalid_argument("input/output sizes must match and output must not alias inputs");
    for (std::size_t i = 0; i < n; ++i) {
        Byte value = gf_mul(column[i], x[0]);
        // Split the diagonal: no per-element index-selection branch.
        for (std::size_t j = 1; j <= i; ++j) value ^= gf_mul(column[i-j], x[j]);
        for (std::size_t j = i+1; j < n; ++j) value ^= gf_mul(row[j-i], x[j]);
        out[i] = value;
    }
    if (counts) { counts->multiplications += n*n; counts->additions += n*(n-1); }
}
inline void brute_force_circulant_into(const Vector& column, const Vector& x, Vector& out, Counts* counts = nullptr) {
    detail::check_size(column.size());
    const auto n = column.size();
    if (x.size() != n || out.size() != n || &out == &x || &out == &column)
        throw std::invalid_argument("input/output sizes must match and output must not alias inputs");
    for (std::size_t i = 0; i < n; ++i) {
        Byte value = gf_mul(column[i], x[0]);
        // N is a power of two; masking replaces the old integer remainder.
        for (std::size_t j = 1; j < n; ++j) value ^= gf_mul(column[(i-j)&(n-1)], x[j]);
        out[i] = value;
    }
    if (counts) { counts->multiplications += n*n; counts->additions += n*(n-1); }
}
inline Vector brute_force_toeplitz(const Vector& column, const Vector& row, const Vector& x, Counts* counts = nullptr) {
    Vector out(column.size());
    brute_force_toeplitz_into(column, row, x, out, counts);
    return out;
}
inline Vector brute_force_circulant(const Vector& column, const Vector& x, Counts* counts = nullptr) {
    Vector out(column.size());
    brute_force_circulant_into(column, x, out, counts);
    return out;
}
} // namespace dynamic_matrix
