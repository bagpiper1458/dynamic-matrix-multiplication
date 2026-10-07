# Derivation, conventions, and operation counts

## Domain

Let N = 2ⁿ, n ≥ 0. Each coefficient is an element of GF(2⁸) in the polynomial basis with modulus 0x11B. Addition and subtraction are both XOR. The algorithm computes T·x, where T is N × N and x has length N. It is not a general dense matrix–matrix multiplication routine.

Toeplitz input is a first column `c` and first row `r`, with `c[0] == r[0]`:

```
T[i,j] = c[i-j]   if i >= j
         r[j-i]   otherwise
```

Circulant input is its first column `a`, and `C[i,j] = a[(i-j) % N]`. For example, `a=[a0,a1,a2,a3]` represents:

```
[a0 a3 a2 a1]
[a1 a0 a3 a2]
[a2 a1 a0 a3]
[a3 a2 a1 a0]
```

## Three recursive products

Every even-sized Toeplitz matrix splits into blocks:

```
T = [A0 A1]     x = [x0]
    [A2 A0]         [x1]

F = A0(x0 XOR x1)
G = (A0 XOR A1)x1
H = (A0 XOR A2)x0

Tx = [F XOR G]
     [F XOR H]
```

Expanding the upper component cancels two copies of A0·x1; expanding the lower cancels two copies of A0·x0. All three smaller matrices are Toeplitz, so the decomposition recurses down to scalar multiplication. This is the identity in the manuscript's correctness argument, implemented directly rather than through its symbolic formula tables.

During preparation, each Toeplitz matrix is represented by its 2N−1 diagonal coefficients, ordered by row-minus-column offset from −(N−1) through N−1. The plan stores three children per internal node. It reuses the **computed value F** for both halves of the output and retains coefficient combinations across repeated calls. It does not perform global symbolic common-subexpression elimination or merge coincident coefficient subtrees. The `PARAMS` pseudocode from the manuscript is therefore replaced by an explicit prepared recursion tree.

## Unit-row-sum circulant case

A circulant matrix has A1 = A2. The matrices A0 XOR A1 and A0 XOR A2 are circulant of half the size, and their row sums equal the original row sum. If XOR(a) = 1, both inherit that property. Propagate a `unit` flag **only into G and H**, never into F, which is generally Toeplitz.

At a flagged scalar leaf, the coefficient is 1 and the result can be copied. There are 2ⁿ = N such leaves, saving exactly N scalar multiplications. Validation rejects a requested unit shortcut if the original coefficient XOR is not 1. A standard integer sum is not the correct test.

## Cost model

M(0)=1 and M(n)=3M(n−1), so M(n)=3ⁿ. With prepared coefficients, the vector XORs obey A(0)=0 and A(n)=3A(n−1)+3·2ⁿ⁻¹, so A(n)=3(3ⁿ−2ⁿ). The unit case saves N multiplications and retains the same additions. No extra opportunistic optimizations for accidental zero/one coefficients are applied.

For this implementation, coefficient preparation obeys P(0)=0 and P(n)=3P(n−1)+2(2ⁿ−1), giving:

```
P(n) = 3·3ⁿ − 4·2ⁿ + 1
```

This includes both G/H coefficient arrays even when they are equal. It is an implementation-specific preparation count, not a claimed minimum. A unit-sum check costs N−1 further XORs in Python's reduction; C++ starts its accumulator at zero and executes N XORs. These checks, allocations, indexing, and bit operations inside field multiplication are excluded from reported field-operation counts. All are included in end-to-end timings.

Each scalar multiplication counts as one M, regardless of internal shifts/XORs. Each byte XOR outside that scalar kernel counts as one field addition A. Baseline dot products start with their first product, so no artificial `0 XOR first_product` is counted.

The prepared plan uses O(3ⁿ) nodes. Both arithmetic and preparation are O(N^(log₂3)); storage and allocation overhead can dominate small instances. Python creates tuples during evaluation. C++ allocates the N-byte output and a 2N−2-byte scratch buffer once per call, then passes slices through recursion without further heap allocations. Its scratch recurrence is W(N)=N+W(N/2), with W(1)=0. The two N/2-byte regions hold the mixed input and F; G/H write directly into the final output halves. A field-operation improvement alone does not establish an execution-time improvement, nor does it prove stronger cryptographic security.
