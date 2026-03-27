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

def S_d(n, t, d):
    """
    S_d = C(n, t-d) * sum_{k=0}^d y_k
    """
    total = sum(y_k(n, t, k) for k in range(d + 1))
    return Fraction(safe_comb(n, t - d), 1) * total

def rhs_S_t(n, m, t):
    """
    RHS of
    S_t = 1/2 * ( 1/2^m + (-1)^t * sum_{j=0}^m (1/2^j) C(n-j-1, m-j) )
    """
    s = Fraction(0, 1)
    for j in range(m + 1):
        s += Fraction(safe_comb(n - j - 1, m - j), 2 ** j)
    return Fraction(1, 2) * (Fraction(1, 2 ** m) + ((-1) ** t) * s)

def rhs_S_t_minus_1(n, m, t):
    """
    RHS of
    S_{t-1} = n/2 * ( 1/2^m + (-1)^(t-1) * sum_{j=0}^m (1/2^j) C(n-j-2, m-j) )
    """
    s = Fraction(0, 1)
    for j in range(m + 1):
        s += Fraction(safe_comb(n - j - 2, m - j), 2 ** j)
    return Fraction(n, 2) * (Fraction(1, 2 ** m) + ((-1) ** (t - 1)) * s)

def verify_identities(max_n=20, print_all=True):
    """
    Verify the two identities for all n <= max_n, under the assumption t = n - m - 1.

    If print_all=True, print lhs and rhs for every tested parameter.
    If print_all=False, print only mismatches.
    """
    all_ok = True

    for n in range(2, max_n + 1):
        for m in range(0, n):
            t = n - m - 1

            print(f"\n{'='*72}")
            print(f"n = {n}, m = {m}, t = {t}")
            print(f"{'='*72}")

            # Check S_t
            lhs1 = S_d(n, t, t)
            rhs1 = rhs_S_t(n, m, t)
            ok1 = (lhs1 == rhs1)

            if print_all or not ok1:
                print("Checking S_t:")
                print(f"  lhs = {lhs1}")
                print(f"  rhs = {rhs1}")
                print(f"  {'OK' if ok1 else 'ERROR'}")

            if not ok1:
                all_ok = False

            # Check S_{t-1}, only when t >= 1
            if t >= 1:
                lhs2 = S_d(n, t, t - 1)
                rhs2 = rhs_S_t_minus_1(n, m, t)
                ok2 = (lhs2 == rhs2)

                if print_all or not ok2:
                    print("Checking S_{t-1}:")
                    print(f"  lhs = {lhs2}")
                    print(f"  rhs = {rhs2}")
                    print(f"  {'OK' if ok2 else 'ERROR'}")

                if not ok2:
                    all_ok = False

    print(f"\n{'='*72}")
    if all_ok:
        print(f"All tests passed for n <= {max_n}.")
    else:
        print(f"Some tests failed for n <= {max_n}.")
    print(f"{'='*72}")

if __name__ == "__main__":
    verify_identities(max_n=20, print_all=True)
