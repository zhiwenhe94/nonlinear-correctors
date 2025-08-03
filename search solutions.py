import itertools
import math


def search_solutions(n, m, t):
    """
    Search for all pairs of vectors Y0 and Y1 satisfying the given constraints:
      - Y_j for j=0..n are non-negative integers
      - For j < n-t, Y_j is either 0 or C(n, j)
      - For j >= n-t, Y_j is computed via the formula:
          Y_{n-t+d} = 2^{n-m-t} * S_d - sum_{j=0}^{n-t-1} T_{j,d} * Y_j
        and must be integer in [0, C(n, n-t+d)]
      - Partition solutions into A0 (Y0=0) and A1 (Y0=1)
      - Combine one vector from A0 and one from A1 such that
          Y0[j] + Y1[j] = C(n, j) for all j
    Returns a list of (Y0, Y1) solution pairs.
    """
    # Precompute binomial coefficients C(n, j)
    C = [math.comb(n, j) for j in range(n + 1)]

    # Compute y_k for k=0..t
    y = [((-1) ** k) * math.comb(n - t + k - 1, k) for k in range(t + 1)]

    # Compute S_d for d=0..t
    S = [math.comb(n, t - d) * sum(y[:d + 1]) for d in range(t + 1)]

    # Compute T_{j,d} for d=0..t and j=0..n-t-1
    T = []
    for d in range(t + 1):
        row = []
        for j in range(n - t):
            total = 0
            for k in range(d + 1):
                idx = j - d + k
                if 0 <= idx <= n - t:
                    total += math.comb(n - t, idx) * y[k]
            row.append((math.comb(n, t - d) / math.comb(n, j)) * total)
        T.append(row)

    # Partitioned solution sets
    A0, A1 = [], []

    # Iterate over choices for Y_j for j in [0..n-t-1]
    for bits in itertools.product([0, 1], repeat=(n - t)):
        Y = [bits[j] * C[j] for j in range(n - t)]
        valid = True

        # Compute Y_{n-t + d} for d in 0..t
        for d in range(t + 1):
            val = (2 ** (n - m - t) * S[d]
                   - sum(T[d][j] * Y[j] for j in range(n - t)))
            if not val.is_integer():
                valid = False
                break
            val = int(val)
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

    # Combine from A0 and A1 to form solution pairs
    solutions = []
    for Y0 in A0:
        for Y1 in A1:
            if all(Y0[j] + Y1[j] == C[j] for j in range(n + 1)):
                solutions.append((Y0, Y1))

    return solutions


if __name__ == "__main__":
    # Example usage: scan n=3,4,5,6 with m=1, t=n-m-1
    for n in [5, 6, 7, 8, 9]:
        m = 1
        t = n - m - 4
        sols = search_solutions(n, m, t)
        print(f"n={n}, m={m}, t={t}, found {len(sols)} solution pairs")
        for idx, (Y0, Y1) in enumerate(sols, start=1):
            print(f" Pair {idx}:")
            print("  Y0 =", Y0)
            print("  Y1 =", Y1)
