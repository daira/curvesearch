# -*- coding: utf-8 -*-
import sys
from multiprocessing import Pool, cpu_count
from traceback import print_exc
from math import ceil
from itertools import takewhile

# Avoid limitations of builtin range.
def range(x, y=None, step=1):
    if y is None:
        (x, y) = (0, x)
    i = x
    while i < y:
        yield i
        i += step


DEFAULT_TWOADICITY = 32

# y^2 = x^3 + 1 never has amicable pairs. <https://arxiv.org/pdf/0912.1831.pdf>
COEFFICIENT_RANGE = range(2, 1000)

TWIST_SECURITY = 0
REQUIRE_PRIMITIVE = False

# Formulae to allow for r:
#   0: x^4 + 3
#   1: x^4 - 3*x^2 + 3
FORMULAE_FOR_r = (0, 1)


# Let x = 2^j.X.
#
# We want:
#   x^4 - x^2 + 1 ~ L bits, so X ~ L/4 - j bits.
#   x^4 - x^2 + 1 = 1 (mod 2^twoadicity)
#   x^4 + 3       = 1 (mod 3^threeadicity)
#
# So
#        2^{4j}.X^4 + 3  = 1 (mod 3^threeadicity)
#   2.(2^{4j-1}.X^4 + 1) = 0 (mod ")
#      2^{4j-1}.X^4 + 1  = 0 (mod ")
#
# Therefore X^4 = -1 / 2^{4j-1} (mod 3^threeadicity).
# We find all solutions for X (mod 3^threeadicity) and use them in the search.
#
# There is another case x^4 - 3*x^2 + 3 = 1 (mod 3^threeadicity) which can be
# handled in a similar way, but requires use of the quadratic formula to
# solve for X (mod 3^threeadicity).

class BruteForce:
    def __init__(self, L, twoadicity):
        j = (twoadicity+1)//2
        threeadicity = int(ceil(twoadicity*log(2, 3)))
        pow2 = 2^twoadicity
        pow3 = 3^threeadicity

        Xbase = int(2^(L/4.0 - j))
        Xend = Xbase * 128

        X_0 = [int(X_) for X_ in (-1 / Mod(2^(4*j - 1), pow3)).nth_root(4, all=True)]
        assert(len(X_0) > 0)

        X_1 = [int(X_) for X_ in sum([(i / Mod(4^j, pow3)).nth_root(2, all=True) for i in (1, 2)], [])]
        assert(len(X_1) > 0)

        Xoffsets = (X_0, X_1)
        #print(Xbase, Xoffsets)

        self.params = (Xbase, Xend, Xoffsets, j, pow2, pow3, L)

        print("Xbase = %s, 2-adicity = %d, 3-adicity = %d" % (format_int(Xbase, 2), twoadicity, threeadicity))

    def run(self, wid, processes):
        (Xbase, Xend, Xoffsets, j, pow2, pow3, L) = self.params

        # Align Xbase for this worker.
        chunksize = pow3*processes
        Xbase = ((Xbase+chunksize-1) // chunksize)*chunksize + wid*pow3

        for Xchunk in range(Xbase, Xend, chunksize):
            for i in FORMULAE_FOR_r:
                rdesc = ("x^4 + 3", "x^4 - 3*x^2 + 3")[i]
                for Xoffset in Xoffsets[i]:
                    x = (Xchunk + Xoffset) << j

                    #print("i = %d, j = %d, x = %s = %s = %d (mod 3)" % (i, j, format_int(x, 2), format_int(x, 3), x%3))
                    if x % 3 != 1:
                        x = x << 1
                        assert(x % 3 == 1)

                    q = x^4 - x^2 + 1
                    if q < 2^L or q % pow2 != 1:
                        continue
                    #print("q = %s, i = %d, len nope = %r, pow2 nope = %r" % (format_int(q, 2), i, q < 2^L, q % pow2 != 1))

                    r = x^4 + 3
                    if i == 1:
                        r = r - 3*x^2
                    if r % pow3 != 1:
                        continue

                    q_xm1sq = q*(x-1)^2
                    assert(q_xm1sq % 3 == 0)
                    p = x + q_xm1sq//3
                    # p is less likely to be prime than q or r, so check p first.
                    if is_pseudoprime(p) and is_pseudoprime(q):
                        sys.stderr.write('.')
                        sys.stderr.flush()
                        if is_prime(r) and is_prime(p) and is_prime(q):
                            yield (x, p, q, r, rdesc)

        sys.stderr.write('<')
        sys.stderr.flush()

def find_nice_curves(*args):
    (strategy, wid, processes) = args
    for (x, p, q, r, rdesc) in strategy.run(wid, processes):
        sys.stderr.write('@')
        sys.stderr.flush()

        #print("\nx = %s\np = %s\nq = %s\nr = %s\n  = %s" % (format_int(x, 2), format_int(p), format_int(q, 2), format_int(r, 3), rdesc))

        cofactor = ((x-1)^2)//3
        (Ep, bp) = find_curve(p, q * cofactor)
        if bp == None: continue
        (Eq, bq) = find_curve(q, r)
        if bq == None: continue
        (Er, br) = find_curve(r, q)
        if br == None: continue

        sys.stdout.write('*')
        sys.stdout.flush()

        primq = (Mod(bq, q).multiplicative_order() == q-1)
        if REQUIRE_PRIMITIVE and not primq: continue
        primr = (Mod(br, r).multiplicative_order() == r-1)
        if REQUIRE_PRIMITIVE and not primr: continue

        (twsecq, twembedq) = twist_security(q, r)
        if twsecq < TWIST_SECURITY: continue
        (twsecr, twembedr) = twist_security(r, q)
        if twsecr < TWIST_SECURITY: continue

        (secp, embedp) = curve_security(p, q)
        (secq, embedq) = curve_security(q, r)
        (secr, embedr) = curve_security(r, q)

        zetaq = GF(q).zeta(3)
        zetaq = min(zetaq, zetaq^2)
        assert(zetaq^3 == Mod(1, q))

        zetar = GF(r).zeta(3)
        Q = Eq.gens()[0]
        zQ = endo(Eq, zetaq, Q)
        if zQ != int(zetar)*Q:
            zetar = zetar^2
            assert(zQ == int(zetar)*Q)
        assert(zetar^3 == Mod(1, r))

        R = Er.gens()[0]
        assert(endo(Er, zetar, R) == int(zetaq)*R)

        embeddivq = (r-1)/embedq
        embeddivr = (q-1)/embedr
        twembeddivq = (2*q + 1 - r)/twembedq
        twembeddivr = (2*r + 1 - q)/twembedr

        yield (x, p, q, r, rdesc, bp, bq, br, zetaq, zetar, primq, primr, secp, secq, secr, twsecq, twsecr,
               embedp, embeddivq, embeddivr, twembeddivq, twembeddivr)

def endo(E, zeta, P):
   (xp, yp) = P.xy()
   return E(zeta*xp, yp)

def find_curve(size, order):
    set_points = set()
    for b in COEFFICIENT_RANGE:
        E = EllipticCurve(GF(size), [0, b])
        points = E.count_points()
        set_points.add(points)
        if points == order:
            return (E, b)
        if len(set_points) == 6:
            break
    return (None, None)

def find_lowest_prime(p):
    for r in Primes():
        if gcd(p-1, r) == 1:
            return r

pi_12 = (pi/12).numerical_approx()

def curve_security(size, order):
    sys.stdout.write('!')
    sys.stdout.flush()
    suborder = factor(order)[-1][0]
    return (log(pi_12 * suborder, 4), embedding_degree(size, suborder))

def twist_security(size, order):
    return curve_security(size, 2*(size+1) - order)

def embedding_degree(size, suborder):
    sys.stdout.write('#')
    sys.stdout.flush()
    assert(gcd(size, suborder) == 1)
    Z_q = Integers(suborder)
    u = Z_q(size)
    d = suborder-1
    V = factor(d)
    for (v, k) in V:
        while d % v == 0:
            if u^(d/v) != 1: break
            d /= v

    return d


def format_int(n, b=None):
    base = b or 10
    if n is None: return 'None'
    n = int(n)
    if n == 0: return '0'
    neg = " " if n > 0 else "-"
    n = abs(n)
    bitlen = n.bit_length()
    nums = []
    while n > 0:
        n, r = divmod(n, base)
        nums.append(str(r))

    adicity = ''
    if b is not None and nums[0] == '1':
        adicity = " (%d-adicity %d)" % (b, 1 + len(list(takewhile(lambda x: x == '0', nums[1:]))))
    return "(%3d bits) %s%s_%d%s" % (bitlen, neg, ''.join(reversed(nums)), base, adicity)


def main():
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    args = [arg for arg in args if not arg.startswith("--")]

    if len(args) < 1:
        print("Usage: sage bls12halogen.sage [--sequential] <min-bitlength> [<min-2adicity>h]\n")
        return

    L          = int(args[0])
    twoadicity = int(args[1]) if len(args) > 1 else DEFAULT_TWOADICITY

    print("Using %d processes." % (processes,))
    pool = Pool(processes=processes)

    strategy = BruteForce(L, twoadicity)

    try:
        for wid in range(processes):
            pool.apply_async(worker, (strategy, wid, processes))

        while True:
            sleep(1000)
    except (KeyboardInterrupt, SystemExit):
        pass
    finally:
        pool.terminate()

def worker(*args):
    try:
        real_worker(*args)
    except (KeyboardInterrupt, SystemExit):
        pass
    except:
        print_exc()

def real_worker(*args):
    for (x, p, q, r, rdesc, bp, bq, br, zetaq, zetar, primq, primr, secp, secq, secr, twsecq, twsecr,
         embedp, embeddivq, embeddivr, twembeddivq, twembeddivr) in find_nice_curves(*args):
        output  = "\n"
        output += "x   = %s\n" % format_int(x, 2)
        output += "p   = %s\n" % format_int(p)
        output += "q   = %s\n" % format_int(q, 2)
        output += "r   = %s\n" % format_int(r, 3)
        output += "    = %s\n" % rdesc
        output += "ζ_q = %s (mod q)\n" % format_int(zetaq)
        output += "ζ_r = %s (mod r)\n" % format_int(zetar)

        output += "Ep/Fp : y^2 = x^3 + %r\n" % (bp,)
        output += "Eq/Fq : y^2 = x^3 + %r\n" % (bq,)
        output += "Er/Fr : y^2 = x^3 + %r\n" % (br,)

        output += "gcd(q-1, %d) = 1\n" % find_lowest_prime(q)
        output += "gcd(r-1, %d) = 1\n" % find_lowest_prime(r)

        output += "%d is %ssquare and %sprimitive in Fq\n" % (bq, "" if Mod(bq, q).is_square() else "non", "" if primq else "non")
        output += "%d is %ssquare and %sprimitive in Fr\n" % (br, "" if Mod(br, r).is_square() else "non", "" if primr else "non")

        output += "Ep Pollard security = %.1f, embedding degree = %d\n" % (secp, embedp)
        output += "Eq Pollard security = %.1f, embedding degree = (r-1)/%d\n" % (secq, embeddivq)
        output += "Er Pollard security = %.1f, embedding degree = (q-1)/%d\n" % (secr, embeddivr)

        output += "Eq twist Pollard security = %.1f, embedding degree = (2q + 1 - r)/%d\n" % (twsecq, twembeddivq)
        output += "Er twist Pollard security = %.1f, embedding degree = (2r + 1 - q)/%d\n" % (twsecr, twembeddivr)

        print(output)  # one syscall to minimize tearing

main()
