#pragma once

#include "dynamic_matrix.hpp"
#include <algorithm>
#include <limits>
#include <unordered_map>

namespace dynamic_matrix {

class IterativeLayout;

// Owns matrix-dependent coefficients; evaluation needs no layout or recursion.
class IterativePlan {
    std::size_t size_, preparation_additions_;
    Vector coefficients_;
    std::vector<std::size_t> identities_;
    friend class IterativeLayout;
    IterativePlan(std::size_t n, std::size_t additions, Vector coefficients,
                  std::vector<std::size_t> identities)
        : size_(n), preparation_additions_(additions),
          coefficients_(std::move(coefficients)), identities_(std::move(identities)) {}
public:
    std::size_t size() const { return size_; }
    std::size_t preparation_additions() const { return preparation_additions_; }
    const Vector& coefficients() const { return coefficients_; }
    // Retained vector capacities only: excludes object, allocator overhead, and RSS.
    std::size_t retained_payload_bytes() const {
        return coefficients_.capacity() + identities_.capacity()*sizeof(std::size_t);
    }
    std::size_t application_buffer_bytes() const { return coefficients_.size() + size_; }

    Vector matvec(const Vector& x, Counts* counts = nullptr) const {
        if (x.size() != size_) throw std::invalid_argument("vector length must match matrix size");
        Vector work(coefficients_.size());
        std::copy(x.begin(), x.end(), work.begin());
        // Expand backwards so writes cannot destroy still-needed binary inputs.
        for (std::size_t stride = 1, groups = size_/2; groups; stride *= 3, groups /= 2) {
            for (std::size_t group = groups; group-- > 0;) {
                const auto source = group*2*stride, target = group*3*stride;
                for (std::size_t i = stride; i-- > 0;) {
                    const Byte a = work[source+i], b = work[source+stride+i];
                    work[target+i] = a ^ b;
                    work[target+stride+i] = b;
                    work[target+2*stride+i] = a;
                }
            }
        }
        std::size_t identity = 0;
        for (std::size_t i = 0; i < work.size(); ++i) {
            if (identity < identities_.size() && i == identities_[identity]) ++identity;
            else work[i] = gf_mul(coefficients_[i], work[i]);
        }
        // Contract forwards, loading all three operands before either write.
        for (std::size_t stride = 1, groups = work.size()/3; groups; stride *= 2, groups /= 3) {
            for (std::size_t group = 0; group < groups; ++group) {
                const auto source = group*3*stride, target = group*2*stride;
                for (std::size_t i = 0; i < stride; ++i) {
                    const Byte a = work[source+i], b = work[source+stride+i], c = work[source+2*stride+i];
                    work[target+i] = a ^ b;
                    work[target+stride+i] = a ^ c;
                }
            }
        }
        if (counts) {
            counts->multiplications += work.size()-identities_.size();
            counts->additions += 3*(work.size()-size_);
            counts->identity_copies += identities_.size();
        }
        return Vector(work.begin(), work.begin()+size_);
    }
};

// Dimension/structure-only coefficient DAG. Reuse across matrices of this shape.
// Canonical operand pairs share symbolic expressions, never coincident byte values.
class IterativeLayout {
    using Pair = std::pair<std::size_t, std::size_t>;
    struct PairHash {
        std::size_t operator()(const Pair& p) const noexcept {
            return p.first ^ (p.second + std::size_t(0x9e3779b9U) + (p.first << 6) + (p.first >> 2));
        }
    };
    std::size_t size_, input_count_;
    bool circulant_;
    std::vector<Pair> operations_;
    std::vector<std::size_t> coefficient_ids_, unit_leaves_;
public:
    explicit IterativeLayout(std::size_t size, bool circulant = false)
        : size_(size), input_count_(0), circulant_(circulant) {
        detail::check_size(size);
        // Bound all ID, buffer-size, and operation-count arithmetic before allocation.
        const auto limit = std::numeric_limits<std::size_t>::max()/16;
        std::size_t leaves = 1;
        for (auto width = size; width > 1; width /= 2) {
            if (leaves > limit/3) throw std::length_error("iterative layout size overflows");
            leaves *= 3;
        }
        input_count_ = circulant ? size : 2*size-1;
        const auto zero = input_count_;
        std::unordered_map<Pair, std::size_t, PairHash> memo;
        auto intern = [&](std::size_t a, std::size_t b) {
            if (a == b) return zero;
            if (a == zero) return b;
            if (b == zero) return a;
            const Pair pair{std::min(a,b), std::max(a,b)};
            const auto found = memo.find(pair);
            if (found != memo.end()) return found->second;
            const auto id = input_count_+1+operations_.size();
            operations_.push_back(pair);
            memo.emplace(pair, id);
            return id;
        };
        std::vector<std::size_t> current(2*size-1);
        for (std::size_t i = 0; i < current.size(); ++i)
            current[i] = circulant ? (i+1)%size : i;
        for (std::size_t width = size, groups = 1; width > 1; width /= 2, groups *= 3) {
            const auto half = width/2, child_width = width-1;
            std::vector<std::size_t> next;
            next.reserve(groups*3*child_width);
            for (std::size_t group = 0; group < groups; ++group) {
                const auto start = group*(2*width-1);
                for (std::size_t i = 0; i < child_width; ++i) next.push_back(current[start+half+i]);
                for (std::size_t i = 0; i < child_width; ++i)
                    next.push_back(intern(current[start+half+i], current[start+i]));
                for (std::size_t i = 0; i < child_width; ++i)
                    next.push_back(intern(current[start+half+i], current[start+width+i]));
            }
            current = std::move(next);
        }
        coefficient_ids_ = std::move(current);
        unit_leaves_.reserve(size);
        for (std::size_t index = 0; index < leaves; ++index) {
            auto code = index;
            bool selected = true;
            for (auto width = size; width > 1; width /= 2, code /= 3) {
                if (code%3 == 0) { selected = false; break; }
            }
            if (selected) unit_leaves_.push_back(index);
        }
    }
    std::size_t size() const { return size_; }
    std::size_t preparation_additions() const { return operations_.size(); }
    std::size_t retained_payload_bytes() const {
        return operations_.capacity()*sizeof(Pair)
             + (coefficient_ids_.capacity()+unit_leaves_.capacity())*sizeof(std::size_t);
    }
    IterativePlan prepare(const Vector& column, const Vector& row = {}, bool unit_row_sum = false) const {
        if (column.size() != size_) throw std::invalid_argument("matrix size must match layout");
        Vector params;
        params.reserve(input_count_+1+operations_.size());
        if (circulant_) {
            if (!row.empty()) throw std::invalid_argument("circulant layouts take a column only");
            if (unit_row_sum) {
                Byte sum = 0;
                for (Byte a : column) sum ^= a;
                if (sum != 1) throw std::invalid_argument("column XOR must equal one");
            }
            params.insert(params.end(), column.begin(), column.end());
        } else {
            if (unit_row_sum) throw std::invalid_argument("unit shortcut requires a circulant layout");
            detail::check_toeplitz(column, row);
            params.insert(params.end(), row.rbegin(), row.rend()-1);
            params.insert(params.end(), column.begin(), column.end());
        }
        params.push_back(0);
        for (const auto& op : operations_) params.push_back(params[op.first] ^ params[op.second]);
        Vector coefficients(coefficient_ids_.size());
        for (std::size_t i = 0; i < coefficients.size(); ++i) coefficients[i] = params[coefficient_ids_[i]];
        std::vector<std::size_t> identities;
        if (unit_row_sum) {
            identities = unit_leaves_;
            for (auto i : identities)
                if (coefficients[i] != 1) throw std::logic_error("unit-row-sum invariant violated");
        }
        return IterativePlan(size_, operations_.size(), std::move(coefficients), std::move(identities));
    }
};
} // namespace dynamic_matrix
