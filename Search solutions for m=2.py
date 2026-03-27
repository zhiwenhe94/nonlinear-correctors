import itertools
from math import comb

def compute_S_T(n, m=2):
    """
    Compute the constants S_d and T_{j,d} for given n and m=2.
    Returns:
      S: list of length t+1, where t = n-m-1
      T: dict mapping (j,d) -> T_{j,d}
    """
    t = n - m - 1
    # 1) compute y_k
    y = [((-1)**k) * comb(n - t + k - 1, k) for k in range(t+1)]
    # 2) compute S_d = C(n, t-d) * sum_{k=0..d} y_k
    S = [comb(n, t-d) * sum(y[:d+1]) for d in range(t+1)]
    # 3) compute T_{j,d}
    T = {}
    for j in range(n - t):          # j = 0..n-t-1
        for d in range(t+1):        # d = 0..t
            numerator = comb(n, t-d)
            denominator = comb(n, j)
            s = 0
            for k in range(d+1):
                idx = j - d + k
                if 0 <= idx <= n - t:
                    s += comb(n - t, idx) * y[k]
            # integer division guaranteed
            T[(j, d)] = numerator * s // denominator
    return S, T

def find_solutions(n, m=2):
    """
    Find all 4-row solutions [Y_0...Y_n] for given n, m=2 and t=n-m-1.
    Algorithm:
      1. Compute S, T.
      2. Enumerate all "single-row" vectors of length n+1 that satisfy
         Y_{n-t+d} = 2^(n-m-t)*S[d] - sum_{j=0..n-t-1} T[j,d]*Y_j,
         with 0 <= Y_j <= C(n,j).
      3. Split single rows into A0 (Y0=0) and A1 (Y0=1).
      4. Choose 3 rows from A0 (with repetition) and 1 from A1,
         check column-sum constraint sum_i Y_{i,j} == C(n,j).
    Returns:
      rows: list of all valid single-row tuples
      solutions: list of lists of 4 tuples (the valid 4-row groups)
    """
    t = n - m - 1
    S, T = compute_S_T(n, m)

    # 2) enumerate single-row solutions
    rows = []
    for y0 in range(comb(n, 0)+1):
        for y1 in range(comb(n, 1)+1):
            for y2 in range(comb(n, 2)+1):
                vals = [y0, y1, y2]
                sol = vals.copy()
                ok = True
                for d in range(t+1):
                    val = (2**(n-m-t) * S[d]
                           - sum(T[(j,d)] * vals[j] for j in range(n-t)))
                    if not (isinstance(val, int) and 0 <= val <= comb(n, t+d)):
                        ok = False
                        break
                    sol.append(val)
                if ok:
                    rows.append(tuple(sol))
    rows = sorted(set(rows))

    # 3) split by Y0
    A0 = [r for r in rows if r[0] == 0]
    A1 = [r for r in rows if r[0] == 1]

    # 4) choose 3 from A0, 1 from A1
    solutions = []
    for combo0 in itertools.combinations_with_replacement(A0, 3):
        for r1 in A1:
            group = list(combo0) + [r1]
            # check column sums
            if all(sum(row[j] for row in group) == comb(n, j)
                   for j in range(n+1)):
                solutions.append(group)

    return rows, solutions

if __name__ == '__main__':
    for n in [3, 4, 5, 6]:
        rows, sols = find_solutions(n)
        print(f'\nn = {n}')
        print('Single-row solutions:')
        for r in rows:
            print(r)
        print(f'\n4-row groups (total {len(sols)}):')
        for sol in sols:
            print(sol)
