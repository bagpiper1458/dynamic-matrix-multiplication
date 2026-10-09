#include "dynamic_matrix.hpp"
#include "iterative_matrix.hpp"
#include "packed_matrix.hpp"
#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>

using namespace dynamic_matrix;

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
Byte oracle(Byte a, Byte b) {
    unsigned product = 0;
    for (int bit = 0; bit < 8; ++bit) if ((b >> bit) & 1) product ^= unsigned(a) << bit;
    for (int degree = 14; degree >= 8; --degree)
        if ((product >> degree) & 1) product ^= 0x11Bu << (degree-8);
    return static_cast<Byte>(product);
}
// Independent small-size oracle: manuscript Algorithm 1/2's set generation.
// SM_B starts at N-1-i and its columns are reversed for the final pairing.
Vector paper_general_formula(const Vector& column, const Vector& row, const Vector& x, bool circular) {
    const auto size = column.size();
    Vector output(size);
    for (std::size_t i = 0; i < size; ++i) {
        std::vector<std::set<int>> a(size), b(size);
        a[0].insert(0);
        b[0].insert(static_cast<int>(size-1-i));
        for (std::size_t j = 1; j < size; ++j) {
            std::size_t step = 1;
            while (2*step <= j) step *= 2;
            const int shift = (i & step) ? static_cast<int>(step) : -static_cast<int>(step);
            a[j] = a[j-step];
            b[j] = b[j-step];
            for (auto offset : a[j-step]) a[j].insert(offset+shift);
            for (auto index : b[j-step]) b[j].insert(index+shift);
        }
        for (std::size_t j = 0; j < size; ++j) {
            Byte coefficient = 0, data = 0;
            for (auto offset : a[j]) {
                const auto wrapped = (offset+static_cast<int>(size))%static_cast<int>(size);
                coefficient ^= circular ? column[wrapped] : offset >= 0 ? column[offset] : row[-offset];
            }
            for (auto index : b[size-1-j]) data ^= x.at(static_cast<std::size_t>(index));
            output[i] ^= oracle(coefficient, data);
        }
    }
    return output;
}
Vector read_vector(std::istream& input, std::size_t n) {
    Vector result(n);
    for (auto& value : result) {
        int parsed = -1;
        if (!(input >> parsed) || parsed < 0 || parsed > 255) throw std::runtime_error("invalid fixture byte");
        value = static_cast<Byte>(parsed);
    }
    return result;
}
template<class F> void rejects(F function) {
    bool caught = false;
    try { function(); } catch (const std::invalid_argument&) { caught = true; }
    require(caught, "invalid input not rejected");
}
int main(int argc, char** argv) {
    try {
        for (int a = 0; a < 256; ++a) for (int b = 0; b < 256; ++b)
            require(gf_mul(static_cast<Byte>(a), static_cast<Byte>(b)) == oracle(static_cast<Byte>(a), static_cast<Byte>(b)), "field mismatch");
        require(prepare_circulant({2,1,1,3}, true).matvec({0xdb,0x13,0x53,0x45}) == Vector({0x8e,0x4d,0xa1,0xbc}), "known vector mismatch");
        std::mt19937 rng(20261007);
        std::size_t power3 = 1;
        for (std::size_t n = 0; n <= 10; ++n, power3 *= 3) {
            const std::size_t size = std::size_t(1) << n;
            Vector column(size), row(size), x(size);
            for (auto& a : column) a = static_cast<Byte>(rng() & 255);
            for (auto& a : row) a = static_cast<Byte>(rng() & 255);
            for (auto& a : x) a = static_cast<Byte>(rng() & 255);
            for (int mode = 0; mode < 3; ++mode) {
                if (mode == 2) {
                    Byte sum = 1;
                    for (std::size_t i = 0; i+1 < size; ++i) sum ^= column[i];
                    column.back() = sum;
                }
                row[0] = column[0];
                auto plan = mode == 0 ? prepare_toeplitz(column, row) : prepare_circulant(column, mode == 2);
                Counts count;
                auto result = plan.matvec(x, &count);
                auto expected = mode == 0 ? toeplitz_matvec(column, row, x) : circulant_matvec(column, x);
                if (size <= 16)
                    require(paper_general_formula(column, row, x, mode != 0) == expected, "literal paper formula mismatch");
                require(result == expected, "random matrix mismatch");
                require(count.multiplications == power3 - (mode == 2 ? size : 0), "multiplication count mismatch");
                require(count.additions == 3*(power3-size), "addition count mismatch");
                require(count.identity_copies == (mode == 2 ? size : 0), "copy count mismatch");
                require(plan.preparation_additions() == 3*power3-4*size+1, "preparation count mismatch");
                require(plan.matvec(Vector(size)) == Vector(size), "zero vector mismatch");
                require(plan.matvec(x) == result, "plan reuse mismatch");
                const IterativeLayout layout(size, mode != 0);
                auto iterative = layout.prepare(column, mode == 0 ? row : Vector{}, mode == 2);
                Counts iterative_count;
                require(iterative.matvec(x, &iterative_count) == expected, "iterative matrix mismatch");
                require(iterative_count.multiplications == count.multiplications &&
                        iterative_count.additions == count.additions &&
                        iterative_count.identity_copies == count.identity_copies, "iterative count mismatch");
                require(iterative.preparation_additions() <= plan.preparation_additions(), "coefficient sharing regressed");
                require(iterative.matvec(Vector(size)) == Vector(size), "iterative zero mismatch");
                require(iterative.matvec(x) == result, "iterative reuse mismatch");
                // Reuse the same dimension schedule with a different matrix.
                const Vector identity_column = [&] { Vector v(size); v[0] = 1; return v; }();
                require(layout.prepare(identity_column, mode == 0 ? identity_column : Vector{}, mode == 2).matvec(x) == x,
                        "layout reuse mismatch");
                for (auto preparation : {Preparation::recursive, Preparation::iterative}) {
                    const auto packed = mode == 0 ? PackedPlan::toeplitz(column, row, preparation) :
                                                   PackedPlan::circulant(column, mode == 2, preparation);
                    require(packed.coefficients() == iterative.coefficients(), "packed preparation coefficient mismatch");
                    require(packed.preparation_additions() == plan.preparation_additions(), "packed preparation count mismatch");
                    Counts rc, ic;
                    require(packed.recursive(x, &rc) == expected, "packed recursive mismatch");
                    require(packed.iterative(x, &ic) == expected, "packed iterative mismatch");
                    require(rc.multiplications == count.multiplications && ic.multiplications == count.multiplications &&
                            rc.additions == count.additions && ic.additions == count.additions &&
                            rc.identity_copies == count.identity_copies && ic.identity_copies == count.identity_copies,
                            "packed operation counts mismatch");
                    Vector out(size), recursive_work(packed.recursive_workspace_bytes()+8, 0xa5),
                           iterative_work(packed.iterative_workspace_bytes()+8, 0xa5), next_x = x;
                    for (auto& value : next_x) value ^= 0x5b;
                    for (const auto& values : {x, Vector(size), next_x}) {
                        const auto truth = mode == 0 ? toeplitz_matvec(column, row, values) : circulant_matvec(column, values);
                        packed.recursive_into(values, out, recursive_work);
                        require(out == truth, "packed reused recursive buffer mismatch");
                        packed.iterative_into(values, out, iterative_work);
                        require(out == truth, "packed reused iterative buffer mismatch");
                        require(std::all_of(recursive_work.end()-8, recursive_work.end(), [](Byte a){return a == 0xa5;}), "recursive guard corrupted");
                        require(std::all_of(iterative_work.end()-8, iterative_work.end(), [](Byte a){return a == 0xa5;}), "iterative guard corrupted");
                    }
                }
                Counts brute_count;
                require((mode == 0 ? brute_force_toeplitz(column, row, x, &brute_count) :
                                     brute_force_circulant(column, x, &brute_count)) == expected, "optimized brute force mismatch");
                require(brute_count.multiplications == size*size && brute_count.additions == size*(size-1), "brute force count mismatch");
            }
        }
        rejects([] { prepare_circulant({}); });
        rejects([] { prepare_circulant({1,2,3}); });
        rejects([] { prepare_circulant({1,1}, true); });
        rejects([] { prepare_toeplitz({1,2}, {2,3}); });
        rejects([] { prepare_toeplitz({1}, {}); });
        rejects([] { prepare_circulant({1,0}).matvec({1}); });
        rejects([] { IterativeLayout bad(0); });
        rejects([] { IterativeLayout bad(3); });
        rejects([] { IterativeLayout(2).prepare({1,2}, {2,3}); });
        rejects([] { IterativeLayout(2).prepare({1,2}); });
        rejects([] { IterativeLayout(2).prepare({1,2}, {1,3}, true); });
        rejects([] { IterativeLayout(2, true).prepare({1}); });
        rejects([] { IterativeLayout(2, true).prepare({1,0}, {1,0}); });
        rejects([] { IterativeLayout(2, true).prepare({1,1}, {}, true); });
        rejects([] { IterativeLayout(2, true).prepare({1,0}).matvec({1}); });
        rejects([] { PackedPlan::circulant({}); });
        rejects([] { PackedPlan::circulant({1,2,3}); });
        rejects([] { PackedPlan::toeplitz({1,2}, {2,3}); });
        rejects([] { PackedPlan::circulant({1,1}, true); });
        rejects([] { PackedPlan::circulant({1,0}).recursive({1}); });
        rejects([] { PackedPlan::circulant({1,0}).iterative({1}); });
        rejects([] { Vector x{1,2}, work(3); PackedPlan::circulant({1,0}).iterative_into(x, x, work); });
        rejects([] { Vector x(4), out(4), work(3); PackedPlan::circulant({1,0,0,0}).recursive_into(x,out,work); });
        rejects([] { Vector x(4), out(4), work(8); PackedPlan::circulant({1,0,0,0}).iterative_into(x,out,work); });
        rejects([] { Vector x{1,2}; brute_force_circulant_into({1,0}, x, x); });
        require(IterativeLayout(4).preparation_additions() == 11, "4x4 Toeplitz PARAMS count");
        require(IterativeLayout(4, true).preparation_additions() == 5, "4x4 circulant PARAMS count");
        require(IterativeLayout(4).prepare({2,7,1,3}, {2,9,11,13}).coefficients() ==
                Vector({2,11,5,9,13,7,3,13,7}), "reviewed 4x4 coefficient ordering mismatch");
        // Plan remains usable after its temporary layout has been destroyed.
        auto independent = IterativeLayout(4, true).prepare({2,1,1,3}, {}, true);
        require(independent.matvec({0xdb,0x13,0x53,0x45}) == Vector({0x8e,0x4d,0xa1,0xbc}), "iterative known vector mismatch");
        std::ifstream fixture(argc > 1 ? argv[1] : "tests/fixtures.txt");
        require(bool(fixture), "cannot open cross-language fixtures");
        std::string line;
        std::size_t cases = 0;
        while (std::getline(fixture, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream input(line);
            char mode;
            std::size_t size;
            require(bool(input >> mode >> size) && size > 0 && size <= 128, "invalid fixture header");
            require(mode == 'T' || mode == 'C' || mode == 'U', "invalid fixture kind");
            auto column = read_vector(input, size);
            auto row = read_vector(input, size);
            auto vector = read_vector(input, size);
            auto expected = read_vector(input, size);
            std::string extra;
            require(!(input >> extra), "extra fixture tokens");
            auto plan = mode == 'T' ? prepare_toeplitz(column, row) : prepare_circulant(column, mode == 'U');
            require(plan.matvec(vector) == expected, "cross-language mismatch");
            auto iterative = IterativeLayout(size, mode != 'T').prepare(column, mode == 'T' ? row : Vector{}, mode == 'U');
            require(iterative.matvec(vector) == expected, "iterative cross-language mismatch");
            for (auto preparation : {Preparation::recursive, Preparation::iterative}) {
                const auto packed = mode == 'T' ? PackedPlan::toeplitz(column, row, preparation) :
                                                PackedPlan::circulant(column, mode == 'U', preparation);
                require(packed.recursive(vector) == expected && packed.iterative(vector) == expected, "packed fixture mismatch");
            }
            require((mode == 'T' ? brute_force_toeplitz(column,row,vector) : brute_force_circulant(column,vector)) == expected,
                    "brute force fixture mismatch");
            ++cases;
        }
        require(cases == 72, "fixture count mismatch");
        std::cout << "PASS: 65536 field products, old/packed recursive/iterative and brute force through N=1024, reused-buffer guards, literal paper formulas through N=16, " << cases << " shared fixtures\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
