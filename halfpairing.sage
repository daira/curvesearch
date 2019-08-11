from math import ceil
from itertools import combinations

# Let E0/Fq : y^2 = x^3 + b0
# Let E1/Fp : y^2 = x^3 + b1

# p and q should each be ~ L bits.

TWOADICITY = 20

def low_hamming_order(l):
    base = 1 << l
    for w in xrange(l-3-TWOADICITY+1):
        for c in combinations(xrange(TWOADICITY, l-3), w):
            yield base + sum([1 << i for i in c])

def find_bn_primes(L):
    # If u = 2^l, then p ~ 36 * 2^4l.
    # Therefore we want l ~ (L - 6)/4.

    l = int(ceil((L - 6)/4.0))
    for u in low_hamming_order(l):
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
    X = format(x, 'b')
    return "0b%s (weight %d)" % (X, sum([int(c) for c in X]))

def find_cycles(L):
    for (p, q, u) in find_bn_primes(L):
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
    print("Usage: sage halfpairing.sage <minbitlength>\n")
else:
    find_cycles(int(sys.argv[1]))
