from fractions import Fraction
import itertools
import math

def search_solutions(n, m, t):
    """
    Search for all pairs of vectors Y0 and Y1 satisfying the given constraints,
    using exact rational arithmetic to avoid floating-point errors.
    """
    # Precompute binomial coefficients C(n, j)
    C = [math.comb(n, j) for j in range(n + 1)]

    # Compute y_k for k=0..t as Fractions
    y = [Fraction(((-1) ** k) * math.comb(n - t + k - 1, k), 1) for k in range(t + 1)]

    # Compute S_d for d=0..t
    S = [Fraction(math.comb(n, t - d), 1) * sum(y[:d + 1]) for d in range(t + 1)]

    # Compute T_{j,d} for d=0..t and j=0..n-t-1 using Fractions
    T = []
    for d in range(t + 1):
        row = []
        for j in range(n - t):
            total = 0
            for k in range(d + 1):
                idx = j - d + k
                if 0 <= idx <= n - t:
                    total += math.comb(n - t, idx) * y[k]
            ratio = Fraction(math.comb(n, t - d), math.comb(n, j))
            row.append(ratio * total)
        T.append(row)

    A0, A1 = [], []

    # Iterate over binary choices for Y_j for j in [0..n-t-1]
    for bits in itertools.product([0, 1], repeat=(n - t)):
        Y = [bits[j] * C[j] for j in range(n - t)]
        valid = True

        # Compute Y_{n-t + d} for d in 0..t
        for d in range(t + 1):
            val = Fraction(2 ** (n - m - t), 1) * S[d] \
                  - sum(T[d][j] * Y[j] for j in range(n - t))
            # Must be an integer
            if val.denominator != 1:
                valid = False
                break
            val = val.numerator
            idx = n - t + d
            if val < 0 or val > C[idx]:
                valid = False
                break
            Y.append(val)

        if not valid or len(Y) != n + 1:
            continue

        # Partition based on Y0
        if Y[0] == 0:
            A0.append(Y.copy())
        elif Y[0] == C[0]:
            A1.append(Y.copy())

    # Combine to form solution pairs
    solutions = []
    for Y0 in A0:
        for Y1 in A1:
            if all(Y0[j] + Y1[j] == C[j] for j in range(n + 1)):
                solutions.append((Y0, Y1))

    return solutions

if __name__ == "__main__":
    # Scan n = 5..20, and all valid t in [0, n-m-1]
    for n in range(5, 21):
        m = 1
        print(f"n={n}, valid t range: 0 to {n-m-1}")
        for t in range(0, n - m):
            sols = search_solutions(n, m, t)
            if sols:
                print(f"  n={n}, t={t}, found {len(sols)} solution pairs")
                # Optionally, print each pair:
                for idx, (Y0, Y1) in enumerate(sols, start=1):
                    print(f"    Pair {idx}:")
                    print(f"      Y0 = {Y0}")
                    print(f"      Y1 = {Y1}")
