# Reconstruction record

## Source

Implemented from the author-provided LaTeX manuscript titled *Improve Computing Efficiency of 2ⁿ × 2ⁿ Toeplitz Matrix and Circulant Matrix Multiplication over Galois Fields by Dynamic Multiplication*. Initially the original sources were unavailable. The author subsequently recovered two closely related C++ formula-generation/validation programs and clarified that PARAMS reuse had only existed as pseudocode. Those programs do not deduplicate the N² product occurrences. The manuscript provides an eight-step scalar multiplication example in C, block identities, operation-count derivations, pseudocode, and historic timing tables. Its bibliography file, figures, compiler settings, benchmark source, and exact hardware configuration were not supplied.

The new implementations follow the block identities and general formula. The [PARAMS derivation](params-derivation.zh-TW.md) defines a concrete interpretation of symbolic keys, block values, preprocessing, and lifetimes. It is new implementation work based on the manuscript and recovered formula code, not recovery of a previously implemented PARAMS optimizer or its performance environment.

## Decisions and discrepancies

1. **Scalar modulus:** use 0x11B, explicitly present in the manuscript's `GFM` function. The separate polynomial x⁸+1 in its discussion concerns cyclic polynomial convolution, not the irreducible scalar-field modulus. For a general N-dimensional cyclic convolution the wraparound polynomial is xᴺ−1, equal to xᴺ+1 in characteristic two. This implementation uses direct matrix indexing and does not depend on those ambiguous general-size statements.
2. **Matrix orientation:** follow the displayed circulant matrix using its first column. Tests include a basis vector that returns that column, guarding against an accidental transpose.
3. **Addition counts:** the manuscript's 2 × 2 derivation counts the vector addition and two result additions, but omits coefficient additions such as a0 XOR a1. Interpret its A count as evaluation with prepared matrix coefficients. Record preparation separately.
4. **Evaluation table:** the table described nearby as the non-unit case lists M=1,5,19,65,...; these equal 3ⁿ−2ⁿ, the unit-row-sum case. Run and label both variants explicitly instead of copying that label.
5. **Introduction typo:** a paragraph repeats 3ⁿ for the unit case, whereas the abstract and later derivation give 3ⁿ−2ⁿ. The executable operation-count tests follow the latter formula and the block proof.
6. **Symbolic algorithms:** the supplied pseudocode does not completely specify initialization or storage behavior. The derived PARAMS specification explicitly chooses aligned dyadic blocks with row-wise partial-sum vectors and symbolic structural keys. Its choices and finite-size symbolic checks are documented; this is not claimed to be a uniquely specified original storage layout.
7. **Timing claims:** original times and percentage improvements are historical manuscript claims, not current measurements. New benchmark outputs include their own metadata. Arithmetic savings do not imply the time-efficiency limit or security claims mentioned in the manuscript.

## Portfolio description

Suggested wording after local verification:

> Reconstructed Python and C++ implementations of a research algorithm for structured matrix–vector multiplication over GF(2⁸). Validated recursive operation counts, built independent arithmetic oracles and cross-language fixtures, and separated preprocessing from evaluation in reproducible benchmarks.

Use a measured speedup in a resume only with the language, compiler settings, matrix size, baseline, and preprocessing policy. The 80% figure from the old resume should not be reused as evidence for this new implementation without a matching new measurement.
