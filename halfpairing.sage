from math import log, floor, ceil
from random import randint
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
    # Therefore we want l ~ (L - lg(36))/4.

    l = int(ceil((L - log(36)/log(2))/4))
    for u in low_hamming_order(l):
        p = 36*u^4 + 36*u^3 + 18*u^2 + 6*u + 1
        if is_pseudoprime(p):
            q = p + 6*u^2
            if is_prime(q) and is_prime(p):
                yield (p, q, u)

def symmetric_range():
    for b in xrange(1, 10000):
        yield b
        yield -b

def find_nice_curve(p, q):
    for b in symmetric_range():
        E = EllipticCurve(GF(p), [0, b])
        if E.count_points() == q:
            return b

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
        print("E0/Fq : y^2 = x^3 + %d" % b0)
        print("E1/Fp : y^2 = x^3 + %d" % b1)
        print("")

find_cycles(400)
