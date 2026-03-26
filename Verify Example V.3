from itertools import product

# ============================================================
# 1. Define the two base Boolean functions on F2^4
# ============================================================

def f1(x):
    a, b, c, d = x
    return (
        (a & b & d)
        ^ (a & b & c)
        ^ (c & d)
        ^ (b & c)
        ^ (a & d)
        ^ (a & c)
        ^ (a & b)
        ^ b ^ c ^ d
    )

def f2(x):
    a, b, c, d = x
    return (
        (a & b & c)
        ^ (a & c & d)
        ^ (a & b)
        ^ (a & d)
        ^ (b & c)
        ^ a ^ c ^ d
    )


# ============================================================
# 2. Basic utility functions
# ============================================================

def dot(u, v):
    """Compute the inner product over F2."""
    return sum(ui * vi for ui, vi in zip(u, v)) % 2

def weight(v):
    """Compute the Hamming weight of a binary vector."""
    return sum(v)

def xor_funcs(func_list):
    """
    Return the XOR-sum of several Boolean functions.
    Each function maps F2^4 to F2.
    """
    def g(x):
        s = 0
        for func in func_list:
            s ^= func(x)
        return s
    return g


# ============================================================
# 3. Walsh transform of a scalar Boolean function on F2^4
# ============================================================

def walsh_transform(func, n=4):
    """
    Compute the Walsh transform of a scalar Boolean function
    func : F2^n -> F2.

    Returns:
        A dictionary W with W[omega] = Walsh coefficient at omega.
    """
    points = list(product([0, 1], repeat=n))
    W = {}

    for omega in points:
        s = 0
        for x in points:
            exponent = (func(x) + dot(omega, x)) % 2
            s += 1 if exponent == 0 else -1
        W[omega] = s

    return W


# ============================================================
# 4. Precompute the four possible block functions:
#    0, f1, f2, f1+f2
# ============================================================

def zero_func(x):
    return 0

base_funcs = {
    (0, 0): zero_func,
    (1, 0): f1,
    (0, 1): f2,
    (1, 1): xor_funcs([f1, f2]),
}

# For each block type, we store:
#   (i) the maximum absolute Walsh coefficient
#   (ii) the Walsh-weight-sum polynomial coefficients
#
# The polynomial is:
#   P_g(z) = sum_{omega in F2^4} W_g(omega) z^{wt(omega)}
#
# If a scalar Boolean function splits as
#   g(x1,...,x5) = g1(x1) + ... + g5(x5),
# then its total Walsh-weight-sum polynomial is
#   P_g(z) = P_{g1}(z) * ... * P_{g5}(z).

base_max_abs = {}
base_poly = {}

for key, func in base_funcs.items():
    W = walsh_transform(func, n=4)

    # Maximum absolute Walsh coefficient
    base_max_abs[key] = max(abs(v) for v in W.values())

    # Group Walsh coefficients by Hamming weight of omega
    coeffs = [0] * 5
    for omega, val in W.items():
        coeffs[weight(omega)] += val
    base_poly[key] = coeffs


# ============================================================
# 5. Polynomial multiplication
# ============================================================

def poly_mul(A, B):
    """
    Multiply two polynomials represented by coefficient lists.
    If
        A(z) = sum_i A[i] z^i,
        B(z) = sum_j B[j] z^j,
    then the result is their convolution.
    """
    C = [0] * (len(A) + len(B) - 1)
    for i, ai in enumerate(A):
        for j, bj in enumerate(B):
            C[i + j] += ai * bj
    return C


# ============================================================
# 6. Describe u·F block by block
# ============================================================

def block_types_for_component(u):
    """
    For a nonzero u in F2^4, compute the block types of u·F.

    Each block type is one of:
        (0,0) -> 0
        (1,0) -> f1
        (0,1) -> f2
        (1,1) -> f1 + f2

    Here F has 5 disjoint blocks x1,...,x5, each in F2^4.
    """
    u1, u2, u3, u4 = u

    # Block x1:
    # u1*f1 + u2*f2 + u3*f1 + u4*f2
    block1 = (u1 ^ u3, u2 ^ u4)

    # Block x2:
    # u1*f1 + u2*f2 + u3*f2 + u4*(f1+f2)
    block2 = (u1 ^ u4, u2 ^ u3 ^ u4)

    # Block x3:
    # u1*f1 + u2*f2 + u3*(f1+f2) + u4*f1
    block3 = (u1 ^ u3 ^ u4, u2 ^ u3)

    # Block x4:
    # u1*f1 + u2*f2
    block4 = (u1, u2)

    # Block x5:
    # u3*f1 + u4*f2
    block5 = (u3, u4)

    return [block1, block2, block3, block4, block5]


# ============================================================
# 7. Compute correction order and vectorial nonlinearity
# ============================================================

def analyze_F():
    """
    Compute:
      - the scalar nonlinearity of each nonzero component u·F,
      - the vectorial nonlinearity N_F,
      - the correction order t of F.

    Since F : F2^20 -> F2^4, the scalar components u·F are Boolean
    functions on 20 variables, but each of them decomposes into 5
    independent 4-variable blocks. This makes the Walsh analysis fast.
    """
    n_total = 20
    m = 4

    nonzero_u = [u for u in product([0, 1], repeat=m) if any(u)]

    detail_N = {}
    weight_sum_data = {}

    for u in nonzero_u:
        block_types = block_types_for_component(u)

        # The Walsh transform factorizes across disjoint blocks:
        # W_{u·F}(omega1,...,omega5) = product_i W_{gi}(omegai)
        #
        # Hence:
        #   max |W_{u·F}| = product_i max |W_{gi}|
        max_abs = 1
        for bt in block_types:
            max_abs *= base_max_abs[bt]

        # Scalar nonlinearity:
        # N = 2^(n-1) - (1/2) * max_abs
        detail_N[u] = 2**(n_total - 1) - max_abs // 2

        # The Walsh-weight-sum polynomial also factorizes:
        # P_{u·F}(z) = product_i P_{gi}(z)
        poly = [1]
        for bt in block_types:
            poly = poly_mul(poly, base_poly[bt])

        weight_sum_data[u] = poly

    # Vectorial nonlinearity
    NF = min(detail_N.values())

    # Correction order:
    # largest t such that for every nonzero u,
    # sum_{wt(omega)=k} W_{u·F}(omega) = 0 for all 0 <= k <= t
    t = -1
    for k in range(n_total + 1):
        ok = all(weight_sum_data[u][k] == 0 for u in nonzero_u)
        if ok:
            t = k
        else:
            break

    return t, NF, detail_N, weight_sum_data


# ============================================================
# 8. Main program
# ============================================================

if __name__ == "__main__":
    t, NF, detail_N, weight_sum_data = analyze_F()

    print("=" * 70)
    print("Analysis of F : F2^20 -> F2^4")
    print("=" * 70)

    print("\nScalar nonlinearities of all nonzero components u·F:")
    for u in sorted(detail_N):
        print(f"u = {u}, N(u·F) = {detail_N[u]}")

    print(f"\nVectorial nonlinearity N_F = {NF}")
    print(f"Correction order t = {t}")

    print("\nWalsh-weight sums for each nonzero component:")
    for u in sorted(weight_sum_data):
        print(f"\nu = {u}")
        for k, val in enumerate(weight_sum_data[u]):
            print(f"  sum_(wt(omega)={k}) W_(u·F)(omega) = {val}")
