import sys
from math import ceil
from itertools import combinations

# Let E0/Fq : y^2 = x^3 + b0
# Let E1/Fp : y^2 = x^3 + b1

# p and q should each be ~ L bits.

DEFAULT_TWOADICITY = 21

def low_hamming_order(l, twoadicity):
    base = 1 << l
    # p-1 is a multiple of 6u, so will have one more trailing zero than u.
    trailing_zeros = twoadicity-1
    for w in xrange(l-trailing_zeros+1):
        for c in combinations(xrange(trailing_zeros, l), w):
            u = base + sum([1 << i for i in c])
            yield u
            yield -u

def find_bn_primes(L, twoadicity):
    # If u = 2^l, then p ~ 36 * 2^4l.
    # Therefore we want l ~ (L - 6)/4.

    l = int(ceil((L - 6)/4.0))
    for u in low_hamming_order(l, twoadicity):
        p = 36*u^4 + 36*u^3 + 18*u^2 + 6*u + 1
        if is_pseudoprime(p):
            q = p + 6*u^2
            if is_prime(q) and is_prime(p):
                yield (p, q, u)

def find_nice_curve(p, q):
    for b in range(1, 10000):
        E = EllipticCurve(GF(p), [0, b])
        if E.count_points() == q:
            return b
    return None

def find_lowest_prime(p):
    for r in Primes():
        if gcd(p-1, r) == 1:
            return r

def format_weight(x):
    X = format(abs(x), 'b')
    return "%s0b%s (weight %d)" % ("-" if x < 0 else "", X, sum([int(c) for c in X]))

def find_cycles(L, twoadicity):
    for (p, q, u) in find_bn_primes(L, twoadicity):
        print("bitlength %d" % len(format(p, 'b')))
        print("p = %s" % format_weight(p))
        print("q = %s" % format_weight(q))
        print("u = %s" % format_weight(u))

        b1 = find_nice_curve(p, q)
        b0 = find_nice_curve(q, p)
        if b0 is None or b1 is None:
            print("No parameters found!")
        else:
            print("E0/Fq : y^2 = x^3 + %d" % b0)
            print("E1/Fp : y^2 = x^3 + %d" % b1)

        print("gcd(p-1, %d) = 1" % find_lowest_prime(p))
        print("gcd(q-1, %d) = 1" % find_lowest_prime(q))
        print("")


if len(sys.argv) <= 1:
    print("Usage: sage halfpairing.sage <min-bitlength> [<min-2adicity>]\n")
else:
    find_cycles(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_TWOADICITY)
