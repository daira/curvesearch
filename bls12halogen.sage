# -*- coding: utf-8 -*-
import sys
from multiprocessing import Pool, cpu_count
from traceback import print_exc
from math import ceil
from itertools import takewhile

if sys.version_info[0] == 2: range = xrange


DEFAULT_TWOADICITY = 32

# y^2 = x^3 + 1 never has amicable pairs. <https://arxiv.org/pdf/0912.1831.pdf>
COEFFICIENT_RANGE = range(2, 1000)

TWIST_SECURITY = 0
REQUIRE_PRIMITIVE = False


# Let x = 2^j.A.X.
#
# We want:
#   x^4 - x^2 + 1 ~ L bits, so A.X ~ L/4 - j bits.
#   x^4 - x^2 + 1 = 1 (mod 2^twoadicity)
#   x^4 + 3       = 1 (mod 3^threeadicity)
#
# So 2^{4j}.A^4.X^4 + 3 = 1 (mod 3^threeadicity)
#    2.(2^{4j-1}.A^4.X^4 + 1) = 0 (mod 3)
#       2^{4j-1}.A^4.X^4 + 1  = 0 (mod 3)
#
# Therefore X^4 = -1 / (2^{4j-1}.A^4) (mod 3^threeadicity).
#
# This has a solution for X with probability around 1/4. If it does, X is
# heuristically about the same length as 3^threeadicity (although it may
# be shorter), and so we have approximately A ~ L/4 - j - twoadicity bits.
#
# There is another case x^4 - 3*x^2 + 3 = 1 (mod 3^adicity) which can be
# handled in a similar way, but requires use of the quadratic formula to
# solve for X.

class BruteForce:
    def __init__(self, L, twoadicity):
        j = (twoadicity+1)//2
        threeadicity = int(ceil(twoadicity*log(2, 3)))

        Alen = (L+3)//4 - j - twoadicity
        Abase = 1 << Alen
        Aend = Abase * 64

        self.pow2 = 2^twoadicity
        self.pow3 = 3^threeadicity
        self.params = (Abase, Aend, j, L)

        print("Alen = %d, 2-adicity = %d, 3-adicity = %d" % (Alen, twoadicity, threeadicity))

    def run(self, wid, processes):
        (Abase, Aend, j, L) = self.params

        # Align Abase for this worker.
        Abase = ((Abase+processes-1) // processes)*processes + wid

        for A in range(Abase, Aend, processes):
            for i in range(2):
                rdesc = ("x^4 + 3", "x^4 - 3*x^2 + 3")[i]
                try:
                    if i == 0:
                        ## We choose k a little smaller than j and filter out values that don't achieve the required 2-adicity.
                        k = j
                        X2 = sqrt(-1 / Mod(2^(4*k - 1) * A^4, self.pow3))
                    else:
                        # Quadratic formula: <http://www.theochem.ru.nl/~pwormer/Knowino/knowino.org/wiki/Quadratic_equation/Advanced.html>
                        # j must be odd
                        k = (j//2)*2 + 1
                        X2 = (3*2^(k-1) + sqrt(Mod(9*2^(2*k - 2) - 2^(4*k + 1), self.pow3))) / 2^(4*k)
                except ZeroDivisionError:
                    continue
                if not X2.is_square():
                    continue
                X = int(sqrt(X2))
                x = (A * X) << k

                # p is less likely to be prime than r, so check p first.
                #print("i = %d, k = %d, A = %s, X = %s, x = %s = %s = %d (mod 3)" % (i, k, format_int(A), format_int(X, 2), format_int(x, 2), format_int(x, 3), x%3))
                if x % 3 == 1:
                    q = x^4 - x^2 + 1
                    if q < 2^L or q % self.pow2 != 1:
                        continue
                    #print("q = %s, A = %d, i = %d, len nope = %r, pow2 nope = %r" % (format_int(q, 2), A, i, q < 2^L, q % self.pow2 != 1))
                    assert(q % self.pow2 == 1)
                    assert((q*(x-1)^2) % 3 == 0)
                    p = x + (q*(x-1)^2)//3
                    r = x^4 + 3
                    if i == 1:
                        r = r - 3*x^2
                    if r % self.pow3 == 1 and is_pseudoprime(p):
                        sys.stderr.write('.')
                        sys.stderr.flush()
                        if is_pseudoprime(q):
                            sys.stderr.write('!')
                            sys.stderr.flush()
                            if is_pseudoprime(r):
                                yield (p, q, r, rdesc, x)

        sys.stderr.write('<')
        sys.stderr.flush()

def find_nice_curves(*args):
    (strategy, wid, processes) = args
    for (p, q, r, rdesc, x) in strategy.run(wid, processes):
        sys.stderr.write('#')
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
        #print(E, points, len(set_points))
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

        #output += "%d is %ssquare and %sprimitive in Fq\n" % (bq, "" if Mod(bq, q).is_square() else "non", "" if primq else "non")
        #output += "%d is %ssquare and %sprimitive in Fr\n" % (br, "" if Mod(br, r).is_square() else "non", "" if primp else "non")

        output += "Ep Pollard security = %.1f, embedding degree = %d\n" % (secp, embedp)
        output += "Eq Pollard security = %.1f, embedding degree = (r-1)/%d\n" % (secq, embeddivq)
        output += "Er Pollard security = %.1f, embedding degree = (q-1)/%d\n" % (secr, embeddivr)

        output += "Eq twist Pollard security = %.1f, embedding degree = (2q + 1 - r)/%d\n" % (twsecq, twembeddivq)
        output += "Er twist Pollard security = %.1f, embedding degree = (2r + 1 - q)/%d\n" % (twsecr, twembeddivr)

        print(output)  # one syscall to minimize tearing

main()
