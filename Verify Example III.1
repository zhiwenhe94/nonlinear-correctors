from itertools import product

# ==========================================
# 1. Define the Boolean functions f1, f2, F
# ==========================================

def f1(x):
    x1, x2, x3, x4 = x
    return (
        (x1 & x2 & x4)
        ^ (x1 & x2 & x3)
        ^ (x3 & x4)
        ^ (x2 & x3)
        ^ (x1 & x4)
        ^ (x1 & x3)
        ^ (x1 & x2)
        ^ x2 ^ x3 ^ x4
    )

def f2(x):
    x1, x2, x3, x4 = x
    return (
        (x1 & x2 & x3)
        ^ (x1 & x3 & x4)
        ^ (x1 & x2)
        ^ (x1 & x4)
        ^ (x2 & x3)
        ^ x1 ^ x3 ^ x4
    )

def F(x):
    return (f1(x), f2(x))


# ==========================================
# 2. Basic utility functions
# ==========================================

def dot(a, b):
    """Compute the inner product over F2."""
    return sum(i * j for i, j in zip(a, b)) % 2

def weight(v):
    """Compute the Hamming weight of a binary vector."""
    return sum(v)

def component_function(u, Fx):
    """
    Compute the scalar Boolean function u·F(x),
    where u is a nonzero vector in F2^m.
    """
    return sum(ui * fi for ui, fi in zip(u, Fx)) % 2


# ==========================================
# 3. Walsh transform
# ==========================================

def walsh_scalar(func, n):
    """
    Compute the Walsh transform of a scalar Boolean function
    func : {0,1}^n -> {0,1}.

    Returns:
        A dictionary W such that W[omega] = Walsh coefficient at omega.
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

def walsh_component_of_F(u, n):
    """
    Compute the Walsh spectrum of the component function u·F.
    """
    return walsh_scalar(lambda x: component_function(u, F(x)), n)


# ==========================================
# 4. Nonlinearity
# ==========================================

def nonlinearity_scalar_from_walsh(W, n):
    """
    Compute the nonlinearity of a scalar Boolean function
    from its Walsh spectrum.
    """
    max_abs = max(abs(v) for v in W.values())
    return 2**(n - 1) - max_abs // 2

def nonlinearity_vectorial(m, n):
    """
    Compute the vectorial nonlinearity of F:
        N_F = min_{u != 0} N_{u·F}

    Returns:
        NF: the vectorial nonlinearity
        data: a dictionary storing N_{u·F} for each nonzero u
    """
    u_list = [u for u in product([0, 1], repeat=m) if any(u)]
    data = {}
    for u in u_list:
        W = walsh_component_of_F(u, n)
        data[u] = nonlinearity_scalar_from_walsh(W, n)
    return min(data.values()), data


# ==========================================
# 5. Correction order
# ==========================================

def correction_order(m, n):
    """
    Compute the largest t such that for every nonzero u in F2^m,
        sum_{wt(omega)=k} W_{u·F}(omega) = 0
    holds for all k = 0, 1, ..., t.

    Returns:
        t: the correction order
        all_weight_sums: Walsh-weight sums for each nonzero u
    """
    u_list = [u for u in product([0, 1], repeat=m) if any(u)]
    all_weight_sums = {}

    for u in u_list:
        W = walsh_component_of_F(u, n)
        sums_by_weight = {}
        for k in range(n + 1):
            sums_by_weight[k] = sum(val for omega, val in W.items() if weight(omega) == k)
        all_weight_sums[u] = sums_by_weight

    t = -1
    for k in range(n + 1):
        ok = all(all_weight_sums[u][k] == 0 for u in u_list)
        if ok:
            t = k
        else:
            break

    return t, all_weight_sums


# ==========================================
# 6. Main program
# ==========================================

if __name__ == "__main__":
    n = 4
    m = 2

    # Print the truth table of F
    print("Truth table of F(x) = (f1(x), f2(x)):")
    for x in product([0, 1], repeat=n):
        print(f"x = {x}, F(x) = {F(x)}")

    print("\n" + "=" * 50)

    # Compute the nonlinearity
    NF, detail_N = nonlinearity_vectorial(m, n)
    print("Nonlinearity of each nonzero component u·F:")
    for u, val in detail_N.items():
        print(f"u = {u}, N(u·F) = {val}")
    print(f"\nVectorial nonlinearity N_F = {NF}")

    print("\n" + "=" * 50)

    # Compute the correction order
    t, weight_sum_data = correction_order(m, n)
    print("Weight sums of Walsh coefficients for each nonzero u:")
    for u, sums in weight_sum_data.items():
        print(f"u = {u}:")
        for k in range(n + 1):
            print(f"  sum_(wt(omega)={k}) W_(u·F)(omega) = {sums[k]}")
    print(f"\nCorrection order t = {t}")
