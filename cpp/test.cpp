#include "dynamic_matrix.hpp"
#include <fstream>
#include <iostream>
#include <random>
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
        for (std::size_t n = 0; n <= 7; ++n, power3 *= 3) {
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
                require(result == expected, "random matrix mismatch");
                require(count.multiplications == power3 - (mode == 2 ? size : 0), "multiplication count mismatch");
                require(count.additions == 3*(power3-size), "addition count mismatch");
                require(count.identity_copies == (mode == 2 ? size : 0), "copy count mismatch");
                require(plan.preparation_additions() == 3*power3-4*size+1, "preparation count mismatch");
                require(plan.matvec(Vector(size)) == Vector(size), "zero vector mismatch");
                require(plan.matvec(x) == result, "plan reuse mismatch");
            }
        }
        rejects([] { prepare_circulant({}); });
        rejects([] { prepare_circulant({1,2,3}); });
        rejects([] { prepare_circulant({1,1}, true); });
        rejects([] { prepare_toeplitz({1,2}, {2,3}); });
        rejects([] { prepare_toeplitz({1}, {}); });
        rejects([] { prepare_circulant({1,0}).matvec({1}); });
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
            ++cases;
        }
        require(cases == 72, "fixture count mismatch");
        std::cout << "PASS: 65536 field products, matrix/count/input tests, " << cases << " shared fixtures\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
