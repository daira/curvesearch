import sys
from math import ceil
from itertools import combinations, chain

# Let E0/Fq : y^2 = x^3 + b0
# Let E1/Fp : y^2 = x^3 + b1

# p and q should each be ~ L bits.

DEFAULT_TWOADICITY = 21

# It is well known that if g is neither a square nor a cube in Fp, then all
# possible group orders an elliptic curve E : y^2 = x^3 + b can have over Fp
# occur as the order of one of the 6 twists with b \in {1, g, g^2, g^3, g^4, g^5}.

# <https://math.stackexchange.com/questions/127251/when-is-a-not-a-cube-mod-p>:
#   If p = 2 (mod 3) then all elements are cubes.
#   If p = 1 (mod 3) then a is a cube iff a^((p-1)/3) = 1.

# <https://cryptojedi.org/papers/pfcpo.pdf> section 2:
# [...] the order of a curve satisfying the norm equation 3V^2 = 4p - t^2 has one
# of the six forms {p+1 +/- t, p+1 +/- (t +/- 3V)/2} [IEEE Std 1363-2000, section
# A.14.2.3, item 6].
#
# We choose 4p = 3V^2 + t^2, where (V-1)/2 and (t-1)/2 are both multiples of 2^twoadicity.
#
# Then 4p = (3(V-1)^2 + 6(V-1) + 3) + ((t-1)^2 + 2(t-1) + 1)
#         = 3(V-1)^2 + 6(V-1) + (t-1)^2 + 2(t-1) + 4
#       p = 3((V-1)/2)^2 + 3(V-1)/2 + ((t-1)/2)^2 + (t-1)/2 + 1
#
# So p-1 will be a multiple of 2^twoadicity, and so will (p+1-t)-1 = (p-1)-(t-1).
#
# We'd also like both p and q to be 1 (mod 3), so that we have efficient endomorphisms
# on both curves. We explicitly check p = 1 (mod 3), and then if t is chosen to be a
# multiple of 3 then (p-1)-(t-1) will be 1 (mod 3) (but we must still check q since it
# is not necessarily that order).

def low_hamming_order(l, twoadicity):
    Vlen = l//2 + 1
    Vbase = 1 << Vlen
    tlen = l//4
    tbase = 1 << tlen
    trailing_zeros = twoadicity+1
    for w in xrange(tlen-trailing_zeros):
        for Vc in combinations(xrange(trailing_zeros, Vlen), w):
            V = Vbase + sum([1 << i for i in Vc]) + 1
            assert(((V-1)/2) % (1<<twoadicity) == 0)
            for tc in chain(combinations(xrange(trailing_zeros, tlen), w),
                            combinations(xrange(trailing_zeros, tlen), w+1)):
                t = tbase + sum([1 << i for i in tc]) + 1
                assert(((t-1)/2) % (1<<twoadicity) == 0)
                if t % 3 != 1:
                    continue
                p4 = 3*V^2 + t^2
                assert(p4 % 4 == 0)
                p = p4//4
                assert(p % (1<<twoadicity) == 1)
                if p % 3 == 1 and is_prime(p):
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    yield p

def find_nonsquare_noncube(p):
    for g_int in range(2, 100):
        g = Mod(g_int, p)
        if g^((p-1)//3) != 1 and g^((p-1)//2) != 1:
            return g
    return None

def find_nice_curves(L, twoadicity):
    for p in low_hamming_order(L-1, twoadicity):
        g = find_nonsquare_noncube(p)
        if g is None: continue

        for i in range(6):
            b1 = g^i
            E1 = EllipticCurve(GF(p), [0, b1])
            q = E1.count_points()
            if q % (1<<twoadicity) == 1 and q % 3 == 1 and is_prime(q):
                b1 = find_coefficient(p, q)
                if b1 is not None:
                    b0 = find_coefficient(q, p)
                    if b0 is not None:
                        yield (p, q, b1, b0)

def find_coefficient(p, q):
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
    for (p, q, b1, b0) in find_nice_curves(L, twoadicity):
        print("")
        print("bitlength %d" % len(format(p, 'b')))
        print("p = %s" % format_weight(p))
        print("q = %s" % format_weight(q))

        print("E0/Fq : y^2 = x^3 + %d" % b0)
        print("E1/Fp : y^2 = x^3 + %d" % b1)

        print("gcd(p-1, %d) = 1" % find_lowest_prime(p))
        print("gcd(q-1, %d) = 1" % find_lowest_prime(q))


if len(sys.argv) <= 1:
    print("Usage: sage amicable.sage <min-bitlength> [<min-2adicity>]\n")
else:
    find_cycles(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_TWOADICITY)
