from itertools import product

# =========================================================
# 1. Define the two scalar Boolean functions f_M and f_M'
# =========================================================

def f_M(x):
    x1, x2, x3, x4, x5 = x
    return (
        (x3 & x4 & x5)
        ^ (x2 & x4 & x5)
        ^ (x2 & x3 & x5)
        ^ (x1 & x3 & x4)
        ^ (x1 & x2 & x4)
        ^ (x1 & x2 & x3)
        ^ (x4 & x5)
        ^ (x3 & x5)
        ^ (x2 & x5)
        ^ (x1 & x4)
        ^ (x1 & x3)
        ^ (x1 & x2)
        ^ x4 ^ x3 ^ x2 ^ x1 ^ 1
    )

def f_M_prime(x):
    x1, x2, x3, x4, x5 = x
    return (
        (x1 & x2 & x3 & x4)
        ^ (x1 & x2 & x3 & x5)
        ^ (x1 & x2 & x3)
        ^ (x1 & x2 & x5)
        ^ (x1 & x3 & x5)
        ^ (x2 & x3 & x4)
        ^ (x2 & x4 & x5)
        ^ (x3 & x4 & x5)
        ^ (x1 & x2)
        ^ (x1 & x3)
        ^ (x1 & x5)
        ^ (x2 & x4)
        ^ (x3 & x4)
        ^ (x4 & x5)
        ^ x1 ^ x2 ^ x3 ^ x5 ^ 1
    )


# ==========================================
# 2. Basic utility functions
# ==========================================

def dot(a, b):
    """Compute the inner product over F2."""
    return sum(i * j for i, j in zip(a, b)) % 2

def weight(v):
    """Compute the Hamming weight of a binary vector."""
    return sum(v)


# ==========================================
# 3. Walsh transform
# ==========================================

def walsh_transform(func, n):
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


# ==========================================
# 4. Nonlinearity
# ==========================================

def nonlinearity_from_walsh(W, n):
    """
    Compute the nonlinearity of a scalar Boolean function
    from its Walsh spectrum:
        N_f = 2^(n-1) - (1/2) * max_{omega} |W_f(omega)|.
    """
    max_abs = max(abs(v) for v in W.values())
    return 2**(n - 1) - max_abs // 2


# ==========================================
# 5. Correction order
# ==========================================

def correction_order_from_walsh(W, n):
    """
    Compute the largest t such that
        sum_{wt(omega)=k} W_f(omega) = 0
    holds for all k = 0, 1, ..., t.

    Returns:
        t: the correction order
        sums_by_weight: the Walsh sums grouped by Hamming weight
    """
    sums_by_weight = {}

    for k in range(n + 1):
        sums_by_weight[k] = sum(val for omega, val in W.items() if weight(omega) == k)

    t = -1
    for k in range(n + 1):
        if sums_by_weight[k] == 0:
            t = k
        else:
            break

    return t, sums_by_weight


# ==========================================
# 6. Analyze one Boolean function
# ==========================================

def analyze_boolean_function(func, n, name="f"):
    """
    Compute and print:
      - truth table
      - number of zeros and ones
      - nonlinearity
      - correction order
      - Walsh-weight sums
    """
    points = list(product([0, 1], repeat=n))

    # Compute the truth table
    values = {x: func(x) for x in points}

    # Compute the zero set and one set
    zeros = [x for x in points if values[x] == 0]
    ones = [x for x in points if values[x] == 1]

    # Compute the Walsh transform
    W = walsh_transform(func, n)

    # Compute the nonlinearity
    Nf = nonlinearity_from_walsh(W, n)

    # Compute the correction order
    t, sums_by_weight = correction_order_from_walsh(W, n)

    # Print the results
    print("=" * 60)
    print(f"Analysis of {name}")
    print("=" * 60)

    print("\nTruth table:")
    for x in points:
        print(f"x = {x}, {name}(x) = {values[x]}")

    print("\nZero set:")
    for x in zeros:
        print(x)
    print(f"Total number of zeros = {len(zeros)}")

    print("\nOne set:")
    for x in ones:
        print(x)
    print(f"Total number of ones = {len(ones)}")

    print("\nWalsh-weight sums:")
    for k in range(n + 1):
        print(f"sum_(wt(omega)={k}) W_{name}(omega) = {sums_by_weight[k]}")

    print(f"\nNonlinearity of {name} = {Nf}")
    print(f"Correction order of {name} = {t}")
    print()


# ==========================================
# 7. Main program
# ==========================================

if __name__ == "__main__":
    n = 5

    analyze_boolean_function(f_M, n, name="f_M")
    analyze_boolean_function(f_M_prime, n, name="f_M_prime")
