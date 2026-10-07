#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace dynamic_matrix {
using Byte = std::uint8_t;
using Vector = std::vector<Byte>;

inline Byte gf_mul(Byte a, Byte b) {
    Byte result = 0;
    for (int i = 0; i < 8; ++i) {
        if (b & 1) result ^= a;
        a = static_cast<Byte>((static_cast<unsigned>(a) << 1) ^ ((a & 0x80) ? 0x11B : 0));
        b >>= 1;
    }
    return result;
}

struct Counts {
    std::size_t multiplications = 0, additions = 0, identity_copies = 0;
};

namespace detail {
inline void check_size(std::size_t n) {
    if (!n || (n & (n-1))) throw std::invalid_argument("size must be a positive power of two");
}
inline void check_toeplitz(const Vector& column, const Vector& row) {
    check_size(column.size());
    if (column.size() != row.size() || column[0] != row[0])
        throw std::invalid_argument("row and column lengths and first elements must match");
}
struct Node {
    std::size_t size;
    Byte coefficient = 0;
    bool identity = false;
    std::unique_ptr<Node> f, g, h;
    explicit Node(std::size_t n) : size(n) {}
};

inline std::unique_ptr<Node> compile(const Vector& d, std::size_t n, bool unit, std::size_t& additions) {
    auto node = std::make_unique<Node>(n);
    if (n == 1) {
        if (unit && d[0] != 1) throw std::logic_error("unit-row-sum invariant violated");
        node->coefficient = d[0];
        node->identity = unit;
        return node;
    }
    const auto half = n/2;
    Vector a0(n-1), g(n-1), h(n-1);
    // d[n-1 + (row-column)] stores each diagonal. Slice A0 starts at half.
    for (std::size_t i = 0; i < n-1; ++i) {
        a0[i] = d[half+i];
        g[i] = a0[i] ^ d[i];
        h[i] = a0[i] ^ d[n+i];
    }
    additions += 2*(n-1);
    node->f = compile(a0, half, false, additions);
    node->g = compile(g, half, unit, additions);
    node->h = compile(h, half, unit, additions);
    return node;
}

inline void evaluate(const Node& node, const Byte* x, Byte* result, Byte* scratch, Counts* counts) {
    if (node.size == 1) {
        if (node.identity) {
            if (counts) ++counts->identity_copies;
            result[0] = x[0];
            return;
        }
        if (counts) ++counts->multiplications;
        result[0] = gf_mul(node.coefficient, x[0]);
        return;
    }
    const auto half = node.size/2;
    Byte* mixed = scratch;
    Byte* f = scratch + half;
    Byte* child_scratch = scratch + node.size;
    for (std::size_t i = 0; i < half; ++i) mixed[i] = x[i] ^ x[half+i];
    evaluate(*node.f, mixed, f, child_scratch, counts);
    evaluate(*node.g, x+half, result, child_scratch, counts);
    evaluate(*node.h, x, result+half, child_scratch, counts);
    for (std::size_t i = 0; i < half; ++i) {
        result[i] ^= f[i];
        result[half+i] ^= f[i];
    }
    if (counts) counts->additions += 3*half;
}
} // namespace detail

class Plan {
    std::size_t size_, preparation_additions_ = 0;
    std::unique_ptr<detail::Node> root_;
    Plan(Vector diagonals, std::size_t n, bool unit) : size_(n) {
        root_ = detail::compile(diagonals, n, unit, preparation_additions_);
    }
    friend Plan prepare_toeplitz(const Vector&, const Vector&);
    friend Plan prepare_circulant(const Vector&, bool);
public:
    std::size_t size() const { return size_; }
    std::size_t preparation_additions() const { return preparation_additions_; }
    Vector matvec(const Vector& vector, Counts* counts = nullptr) const {
        if (vector.size() != size_) throw std::invalid_argument("vector length must match matrix size");
        Vector result(size_), scratch(2*size_-2);
        detail::evaluate(*root_, vector.data(), result.data(), scratch.data(), counts);
        return result;
    }
};

inline Plan prepare_toeplitz(const Vector& column, const Vector& row) {
    detail::check_toeplitz(column, row);
    Vector diagonals(row.rbegin(), row.rend()-1);
    diagonals.insert(diagonals.end(), column.begin(), column.end());
    return Plan(std::move(diagonals), column.size(), false);
}

inline Plan prepare_circulant(const Vector& column, bool unit_row_sum = false) {
    detail::check_size(column.size());
    Byte sum = 0;
    if (unit_row_sum) {
        for (Byte a : column) sum ^= a;
        if (sum != 1) throw std::invalid_argument("column XOR must equal one");
    }
    Vector diagonals(column.begin()+1, column.end());
    diagonals.insert(diagonals.end(), column.begin(), column.end());
    return Plan(std::move(diagonals), column.size(), unit_row_sum);
}

inline Vector toeplitz_matvec(const Vector& column, const Vector& row, const Vector& x, Counts* counts = nullptr) {
    detail::check_toeplitz(column, row);
    const auto n = column.size();
    if (x.size() != n) throw std::invalid_argument("vector length must match matrix size");
    Vector result(n);
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = gf_mul(column[i], x[0]);
        for (std::size_t j = 1; j < n; ++j)
            result[i] ^= gf_mul(i >= j ? column[i-j] : row[j-i], x[j]);
    }
    if (counts) { counts->multiplications += n*n; counts->additions += n*(n-1); }
    return result;
}

inline Vector circulant_matvec(const Vector& column, const Vector& x, Counts* counts = nullptr) {
    detail::check_size(column.size());
    const auto n = column.size();
    if (x.size() != n) throw std::invalid_argument("vector length must match matrix size");
    Vector result(n);
    for (std::size_t i = 0; i < n; ++i) {
        result[i] = gf_mul(column[i], x[0]);
        for (std::size_t j = 1; j < n; ++j)
            result[i] ^= gf_mul(column[(i+n-j)%n], x[j]);
    }
    if (counts) { counts->multiplications += n*n; counts->additions += n*(n-1); }
    return result;
}
} // namespace dynamic_matrix
