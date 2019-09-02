import sys
from math import ceil
from itertools import combinations

# Let E0/Fq : y^2 = x^3 + b0
# Let E1/Fp : y^2 = x^3 + b1

# p and q should each be ~ L bits.

DEFAULT_TWOADICITY = 21
DEFAULT_STRETCH = 0

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

def low_hamming_order(L, twoadicity):
    Vlen = (L-1)//2 + 1
    Vbase = 1 << Vlen
    tlen = (L-1)//4
    tbase = 1 << tlen
    trailing_zeros = twoadicity+1
    for w in xrange(tlen-trailing_zeros):
        for Vc in combinations(xrange(trailing_zeros, Vlen), w):
            V = Vbase + sum([1 << i for i in Vc]) + 1
            assert(((V-1)/2) % (1<<twoadicity) == 0)
            for tw in xrange(1, w+1):
                for tc in combinations(xrange(trailing_zeros, tlen), tw):
                    t = tbase + sum([1 << i for i in tc]) + 1
                    assert(((t-1)/2) % (1<<twoadicity) == 0)
                    if t % 3 != 1:
                        continue
                    p4 = 3*V^2 + t^2
                    assert(p4 % 4 == 0)
                    p = p4//4
                    assert(p % (1<<twoadicity) == 1)
                    if p % 3 == 1 and is_prime(p):
                        yield p

def near_powerof2_order(L, twoadicity):
    trailing_zeros = twoadicity+1
    Vbase = isqrt((1<<(L+2))//3) >> trailing_zeros
    for Voffset in symmetric_range(10000):
        V = ((Vbase + Voffset) << trailing_zeros) + 1
        assert(((V-1)/2) % (1 << twoadicity) == 0)
        tmp = (1<<(L+2)) - 3*V^2
        if tmp < 0: continue
        tbase = isqrt(tmp) >> trailing_zeros
        for toffset in symmetric_range(10000):
            t = ((tbase + toffset) << trailing_zeros) + 1
            assert(((t-1)/2) % (1<<twoadicity) == 0)
            if t % 3 != 1:
                continue
            p4 = 3*V^2 + t^2
            assert(p4 % 4 == 0)
            p = p4//4
            assert(p % (1<<twoadicity) == 1)
            if p % 3 == 1 and is_prime(p):
                yield p

def find_nonsquare_noncube(p):
    for g_int in xrange(2, 100):
        g = Mod(g_int, p)
        if g^((p-1)//3) != 1 and g^((p-1)//2) != 1:
            return g
    return None

def symmetric_range(n):
    for i in xrange(n):
        yield -i
        yield i+1

def find_nice_curves(strategy, L, twoadicity, stretch):
    for p in strategy(L, max(0, twoadicity-stretch)):
        sys.stdout.write('.')
        sys.stdout.flush()
        if p % (1<<twoadicity) != 1: continue
        gp = find_nonsquare_noncube(p)
        if gp is None: continue

        for i in xrange(6):
            b1 = gp^i
            E1 = EllipticCurve(GF(p), [0, b1])
            q = E1.count_points()
            if q % (1<<twoadicity) == 1 and q % 3 == 1 and is_prime(q):
                b1 = find_coefficient(p, q)
                if b1 is not None:
                    b0 = find_coefficient(q, p)
                    if b0 is not None:
                        gq = find_nonsquare_noncube(q)
                        alpha = gq^((q-1)//3)
                        assert(alpha^3 == Mod(1, q))
                        beta  = gp^((p-1)//3)
                        assert(beta^3  == Mod(1, p))
                        yield (p, q, b1, b0, alpha, beta)

def find_coefficient(p, q):
    for b in xrange(1, 10000):
        E = EllipticCurve(GF(p), [0, b])
        if E.count_points() == q:
            return b
    return None

def find_lowest_prime(p):
    for r in Primes():
        if gcd(p-1, r) == 1:
            return r

def format_weight(x, detail=True):
    X = format(abs(x), 'b')
    if detail:
        assert(X.endswith('1'))
        detailstr = " (weight %d, 2-adicity %d)" % (sum([int(c) for c in X]), len(X) - len(X[:-1].rstrip('0')))
    else:
        detailstr = ""

    return "%s0b%s%s" % ("-" if x < 0 else "", X, detailstr)

def find_and_print(strategy, L, twoadicity, stretch):
    for (p, q, b1, b0, alpha, beta) in find_nice_curves(strategy, L, twoadicity, stretch):
        print("")
        print("bitlength %d" % len(format(p, 'b')))
        print("p = %s" % format_weight(p))
        print("q = %s" % format_weight(q))
        print("alpha = %s (mod q)" % format_weight(int(alpha), detail=False))
        print("beta  = %s (mod p)" % format_weight(int(beta),  detail=False))

        print("E0/Fq : y^2 = x^3 + %d (%ssquare)" % (b0, "" if Mod(b0, q).is_square() else "non"))
        print("E1/Fp : y^2 = x^3 + %d (%ssquare)" % (b1, "" if Mod(b1, p).is_square() else "non"))

        print("gcd(p-1, %d) = 1" % find_lowest_prime(p))
        print("gcd(q-1, %d) = 1" % find_lowest_prime(q))


strategy = near_powerof2_order if "--nearpowerof2" in sys.argv[1:] else low_hamming_order
args = [arg for arg in sys.argv[1:] if not arg.startswith("--")]

if len(args) < 1:
    print("Usage: sage amicable.sage [--nearpowerof2] <min-bitlength> [<min-2adicity> [<stretch]]\n")
else:
    find_and_print(strategy,
                   int(args[0]),
                   int(args[1]) if len(args) > 1 else DEFAULT_TWOADICITY,
                   int(args[2]) if len(args) > 2 else DEFAULT_STRETCH)
