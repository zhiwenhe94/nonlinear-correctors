from itertools import product

# =========================================================
# 1. Define the vectorial Boolean function F = (f1, f2, f3)
# =========================================================

def f1(x):
    x1, x2, x3, x4, x5 = x
    return (
        (x1 & x2 & x3 & x5)
        ^ (x2 & x3 & x4 & x5)
        ^ (x1 & x2 & x3)
        ^ (x2 & x4 & x5)
        ^ (x1 & x4 & x5)
        ^ (x2 & x3 & x4)
        ^ (x2 & x3 & x5)
        ^ (x3 & x4 & x5)
        ^ (x1 & x2 & x5)
        ^ (x1 & x2 & x4)
        ^ (x1 & x2)
        ^ (x1 & x3)
        ^ (x1 & x4)
        ^ (x1 & x5)
        ^ x1 ^ x2 ^ x3 ^ x4 ^ x5
    )

def f2(x):
    x1, x2, x3, x4, x5 = x
    return (
        (x1 & x2 & x3 & x4)
        ^ (x1 & x2 & x4 & x5)
        ^ (x1 & x2 & x3)
        ^ (x3 & x4 & x5)
        ^ (x1 & x2 & x4)
        ^ (x1 & x4 & x5)
        ^ (x3 & x4)
        ^ (x2 & x3)
        ^ (x1 & x5)
        ^ (x3 & x5)
        ^ (x4 & x5)
        ^ (x2 & x4)
        ^ x1 ^ x2 ^ x3 ^ x4
    )

def f3(x):
    x1, x2, x3, x4, x5 = x
    return (
        (x1 & x2 & x3 & x4)
        ^ (x1 & x3 & x4 & x5)
        ^ (x1 & x2 & x3 & x5)
        ^ (x2 & x3 & x4 & x5)
        ^ (x2 & x3 & x5)
        ^ (x1 & x2 & x5)
        ^ (x1 & x2 & x3)
        ^ (x3 & x4 & x5)
        ^ (x1 & x2)
        ^ (x1 & x4)
        ^ x1 ^ x2 ^ x5
    )

def F(x):
    return (f1(x), f2(x), f3(x))


# =========================================================
# 2. Define H(y) = (h1(y), h2(y))
# =========================================================

def h1(y):
    y1, y2, y3 = y
    return y1 ^ y2 ^ (y1 & y2) ^ (y1 & y3)

def h2(y):
    y1, y2, y3 = y
    return 1 ^ y1 ^ y2 ^ (y2 & y3)

def HF(x):
    """
    Compute H(F(x)).
    """
    y = F(x)
    return (h1(y), h2(y))


# =========================================================
# 3. Basic utility functions
# =========================================================

def dot(a, b):
    """
    Compute the inner product over F2.
    """
    return sum(i * j for i, j in zip(a, b)) % 2

def weight(v):
    """
    Compute the Hamming weight of a binary vector.
    """
    return sum(v)

def component_function(u, Gx):
    """
    Compute the scalar Boolean function u·G(x),
    where u is a nonzero vector in F2^m.
    """
    return sum(ui * gi for ui, gi in zip(u, Gx)) % 2


# =========================================================
# 4. Walsh transform
# =========================================================

def walsh_scalar(func, n):
    """
    Compute the Walsh transform of a scalar Boolean function
    func : {0,1}^n -> {0,1}.

    Returns:
        A dictionary W such that W[omega] is the Walsh coefficient
        at omega.
    """
    W = {}
    points = list(product([0, 1], repeat=n))

    for omega in points:
        s = 0
        for x in points:
            exponent = (func(x) + dot(omega, x)) % 2
            s += 1 if exponent == 0 else -1
        W[omega] = s

    return W

def walsh_component_of_vectorial_function(u, G, n):
    """
    Compute the Walsh spectrum of the scalar component u·G.
    """
    return walsh_scalar(lambda x: component_function(u, G(x)), n)


# =========================================================
# 5. Nonlinearity
# =========================================================

def nonlinearity_scalar_from_walsh(W, n):
    """
    Compute the nonlinearity of a scalar Boolean function
    from its Walsh spectrum:
        N_f = 2^(n-1) - (1/2) * max_{omega} |W_f(omega)|.
    """
    max_abs = max(abs(v) for v in W.values())
    return 2**(n - 1) - max_abs // 2


# =========================================================
# 6. Correction order of a vectorial Boolean function
# =========================================================

def correction_order_and_nonlinearity(G, n, m):
    """
    Compute the correction order and vectorial nonlinearity of
    a vectorial Boolean function G : F2^n -> F2^m.

    Returns:
        t: correction order
        NF: vectorial nonlinearity
        detail_N: scalar nonlinearity of each nonzero component u·G
        weight_sum_data: Walsh-weight sums for each nonzero u
    """
    u_list = [u for u in product([0, 1], repeat=m) if any(u)]

    detail_N = {}
    weight_sum_data = {}

    for u in u_list:
        W = walsh_component_of_vectorial_function(u, G, n)

        # Compute scalar nonlinearity of u·G
        detail_N[u] = nonlinearity_scalar_from_walsh(W, n)

        # Group Walsh coefficients by Hamming weight of omega
        sums_by_weight = {}
        for k in range(n + 1):
            sums_by_weight[k] = sum(
                val for omega, val in W.items() if weight(omega) == k
            )
        weight_sum_data[u] = sums_by_weight

    # Vectorial nonlinearity
    NF = min(detail_N.values())

    # Correction order
    t = -1
    for k in range(n + 1):
        ok = all(weight_sum_data[u][k] == 0 for u in u_list)
        if ok:
            t = k
        else:
            break

    return t, NF, detail_N, weight_sum_data


# =========================================================
# 7. Display the results
# =========================================================

def print_analysis(name, G, n, m):
    """
    Print the correction order, vectorial nonlinearity,
    scalar nonlinearities of all nonzero components,
    and Walsh-weight sums.
    """
    t, NF, detail_N, weight_sum_data = correction_order_and_nonlinearity(G, n, m)

    print("=" * 70)
    print(f"Analysis of {name}")
    print("=" * 70)

    print("\nScalar nonlinearities of all nonzero components u·G:")
    for u in sorted(detail_N):
        print(f"u = {u}, N(u·G) = {detail_N[u]}")

    print(f"\nVectorial nonlinearity of {name} = {NF}")
    print(f"Correction order of {name} = {t}")

    print("\nWalsh-weight sums:")
    for u in sorted(weight_sum_data):
        print(f"u = {u}")
        for k in range(n + 1):
            print(f"  sum_(wt(omega)={k}) W_(u·G)(omega) = {weight_sum_data[u][k]}")
        print()


# =========================================================
# 8. Main program
# =========================================================

if __name__ == "__main__":
    n = 5

    # Analyze F : F2^5 -> F2^3
    print_analysis("F", F, n, m=3)

    # Analyze H(F(x)) : F2^5 -> F2^2
    print_analysis("H(F(x))", HF, n, m=2)
