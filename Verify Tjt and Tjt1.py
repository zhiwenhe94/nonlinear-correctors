from fractions import Fraction
import math

def safe_comb(n, k):
    """
    Return C(n, k), with the convention that it is 0
    if k < 0, k > n, or n < 0.
    """
    if n < 0 or k < 0 or k > n:
        return 0
    return math.comb(n, k)

def y_k(n, t, k):
    """
    y_k = (-1)^k * C(n - t + k - 1, k)
    """
    return Fraction(((-1) ** k) * safe_comb(n - t + k - 1, k), 1)

def T_jd(n, t, j, d):
    """
    T_{j,d} = [C(n, t-d) / C(n, j)] * sum_{k=0}^d C(n-t, j-d+k) * y_k
    computed exactly using Fraction.
    """
    numerator_factor = safe_comb(n, t - d)
    denominator_factor = safe_comb(n, j)

    if denominator_factor == 0:
        raise ValueError(f"C({n},{j}) = 0, so T_{{j,d}} is not defined for j={j}.")

    total = Fraction(0, 1)
    for k in range(d + 1):
        total += safe_comb(n - t, j - d + k) * y_k(n, t, k)

    return Fraction(numerator_factor, denominator_factor) * total

def verify_formulas(max_n=20, print_all=True):
    """
    Verify, for t = n - m - 1, the identities
        T_{j,t}   = (-1)^t C(n-j-1, m-j),
        T_{j,t-1} = (-1)^(t-1) (n-j) C(n-j-2, m-j),
    for all 0 <= j <= m.

    Parameters
    ----------
    max_n : int
        Maximum n to test.
    print_all : bool
        If True, print lhs and rhs for every tested parameter.
        If False, print only mismatches.
    """
    all_ok = True

    for n in range(3, max_n + 1):
        for m in range(0, n):
            t = n - m - 1

            print(f"\n{'='*72}")
            print(f"n = {n}, m = {m}, t = {t}")
            print(f"{'='*72}")

            # Check T_{j,t}
            print("\nChecking T_{j,t}:")
            for j in range(0, m + 1):
                lhs1 = T_jd(n, t, j, t)
                rhs1 = Fraction(((-1) ** t) * safe_comb(n - j - 1, m - j), 1)
                ok1 = (lhs1 == rhs1)

                if print_all or not ok1:
                    print(f"  j = {j:2d} | lhs = {lhs1} | rhs = {rhs1} | {'OK' if ok1 else 'ERROR'}")

                if not ok1:
                    all_ok = False

            # Check T_{j,t-1}, only when t >= 1
            if t >= 1:
                print("\nChecking T_{j,t-1}:")
                for j in range(0, m + 1):
                    lhs2 = T_jd(n, t, j, t - 1)
                    rhs2 = Fraction(((-1) ** (t - 1)) * (n - j) * safe_comb(n - j - 2, m - j), 1)
                    ok2 = (lhs2 == rhs2)

                    if print_all or not ok2:
                        print(f"  j = {j:2d} | lhs = {lhs2} | rhs = {rhs2} | {'OK' if ok2 else 'ERROR'}")

                    if not ok2:
                        all_ok = False

    print(f"\n{'='*72}")
    if all_ok:
        print(f"All tests passed for n <= {max_n}.")
    else:
        print(f"Some tests failed for n <= {max_n}.")
    print(f"{'='*72}")

if __name__ == "__main__":
    verify_formulas(max_n=20, print_all=True)
