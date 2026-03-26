import itertools
from math import comb
from functools import lru_cache

def compute_S_T(n, m=3):
    t = n - m - 1
    y = [((-1) ** k) * comb(n - t + k - 1, k) for k in range(t + 1)]
    S = [comb(n, t - d) * sum(y[:d + 1]) for d in range(t + 1)]

    T = {}
    for j in range(n - t):          # j = 0, ..., n-t-1
        for d in range(t + 1):      # d = 0, ..., t
            numerator = comb(n, t - d)
            denominator = comb(n, j)
            s = 0
            for k in range(d + 1):
                idx = j - d + k
                if 0 <= idx <= n - t:
                    s += comb(n - t, idx) * y[k]
            T[(j, d)] = numerator * s // denominator
    return S, T


def enumerate_single_rows(n, m=3):
    t = n - m - 1
    if t < 0:
        return []

    S, T = compute_S_T(n, m)
    free_count = n - t   # m=3 时为 4

    rows = []
    ranges = [range(comb(n, j) + 1) for j in range(free_count)]

    for vals in itertools.product(*ranges):
        sol = list(vals)
        ok = True

        for d in range(t + 1):
            idx = n - t + d
            val = 2 * S[d] - sum(T[(j, d)] * vals[j] for j in range(free_count))
            if not (0 <= val <= comb(n, idx)):
                ok = False
                break
            sol.append(val)

        if ok and len(sol) == n + 1:
            rows.append(tuple(sol))

    return sorted(set(rows))


def col_sum(group, n):
    return tuple(sum(row[j] for row in group) for j in range(n + 1))


def target_vector(n):
    return tuple(comb(n, j) for j in range(n + 1))


def find_two_solutions(n, m=3, limit=2):
    rows = enumerate_single_rows(n, m)
    total_rows = 2 ** m   # m=3 时为 8

    A0 = [r for r in rows if r[0] == 0]
    A1 = [r for r in rows if r[0] == 1]

    solutions = []
    tgt = target_vector(n)

    A0 = sorted(A0, key=sum, reverse=True)

    @lru_cache(None)
    def possible(start, k, target):
        target = tuple(target)

        if k == 0:
            return all(x == 0 for x in target)
        if start >= len(A0):
            return False
        if any(x < 0 for x in target):
            return False

        for i in range(start, len(A0)):
            row = A0[i]
            if all(row[j] <= target[j] for j in range(n + 1)):
                new_target = tuple(target[j] - row[j] for j in range(n + 1))
                if possible(i, k - 1, new_target):
                    return True
        return False

    def dfs(start, k, target, chosen):
        if len(solutions) >= limit:
            return True

        if k == 0:
            if all(x == 0 for x in target):
                group = chosen[:]
                if col_sum(group, n) == tgt:
                    solutions.append(group)
                    return len(solutions) >= limit
            return False

        for i in range(start, len(A0)):
            row = A0[i]
            if all(row[j] <= target[j] for j in range(n + 1)):
                new_target = tuple(target[j] - row[j] for j in range(n + 1))
                if possible(i, k - 1, new_target):
                    chosen.append(row)
                    stop = dfs(i, k - 1, new_target, chosen)
                    chosen.pop()
                    if stop:
                        return True
        return False

    for r1 in A1:
        remaining = tuple(tgt[j] - r1[j] for j in range(n + 1))
        possible.cache_clear()
        dfs(0, total_rows - 1, remaining, [r1])
        if len(solutions) >= limit:
            break

    return rows, A0, A1, solutions


if __name__ == '__main__':
    m = 3
    for n in range(4, 12):
        rows, A0, A1, sols = find_two_solutions(n, m=m, limit=2)

        print(f"\n{'='*70}")
        print(f"n = {n}, m = {m}, t = {n-m-1}")
        print(f"single-row solutions: {len(rows)}")
        print("All single-row solutions:")
        for r in rows:
            print(r)

        print(f"\nA0 size = {len(A0)}")
        print("A0 rows:")
        for r in A0:
            print(r)

        print(f"\nA1 size = {len(A1)}")
        print("A1 rows:")
        for r in A1:
            print(r)

        print(f"\nvalid solutions found: {len(sols)}")
        for idx, sol in enumerate(sols, 1):
            print(f"\nSolution {idx}:")
            for row in sol:
                print(row)
            print("column sum =", col_sum(sol, n))
            print("target     =", target_vector(n))
            print("check      =", col_sum(sol, n) == target_vector(n))
