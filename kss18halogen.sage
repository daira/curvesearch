# -*- coding: utf-8 -*-
import sys
from multiprocessing import Pool, cpu_count
from traceback import print_exc
from math import ceil
from itertools import takewhile, starmap

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

lg343 = log(343, 2).numerical_approx()

def solveSextic(a, b, c, n):
    # return solutions for a*x^6 + b*x^3 + c = 0 (mod pow3)
    d = Mod(b, n)^2 - 4*a*c
    if not d.is_square(): return []
    e = sqrt(d)
    sols2 = [(-b + e)/(2*a), (-b - e)/(2*a)]
    sols6 = [int(X_) for X_ in sum([s.nth_root(3, all=True) for s in sols2], [])]
    for x_ in sols6:
        assert((a*x_^6 + b*x_^3 + c) % n == 0)
    return sols6


class BruteForce:
    def __init__(self, L, twoadicity):
        j = (twoadicity+2)//3
        threeadicity = int(ceil(twoadicity*log(2, 3)))
        pow2 = 2^twoadicity
        pow3 = 3^threeadicity

        # Let x = 2^j.X.
        #
        # We want:
        #   343*x^6 + 37*x^3 + 1 ~ L bits, so X ~ (L - lg(343))/6 - j bits.
        #   343*x^6 + 37*x^3 + 1 = 1 (mod 2^twoadicity)
        #   343*x^6 +  u*x^3 + v = 1 (mod 3^threeadicity)

        Xbase = int(2^((L - lg343)/6.0 - j))
        Xend = Xbase * 128

        solutions = [
            (54, 3, #    343*x^6          + 54*x^3          + 3  = q + 1 + T
                    #                                            = 1 (mod 3^threeadicity)
                    # 2*(343*2^{6j-1}*X^6 + 54*2^{3j-1}*X^3 + 1) = 0 (mod ")
                    #    343*2^{6j-1}*X^6 + 54*2^{3j-1}*X^3 + 1  = 0 (mod ")
                    #
                    solveSextic(343*2^(6*j-1), 54*2^(3*j-1), 1, pow3)
            ),
            (20, 1, #    343*x^6          + 20*x^3          + 1  = q + 1 - T
                    #                                            = 1 (mod 3^threeadicity)
                    #    343*2^{6j}*X^6   + 20*2^{3j}*X^3   + 1  = 1 (mod ")
                    #   (343*2^{6j}*X^3   + 20*2^{3j})*X^3       = 0 (mod ")
                    #
                    # X = 0 or X^3 = -20 / 343.2^{3j} (mod 3^threeadicity)
                    [0] + [int(X_) for X_ in (-20 / (343*Mod(2^(3*j), pow3))).nth_root(3, all=True)]
            ),
            (74, 4, #    343*x^6          + 74*x^3          + 4  = q + 1 + (T+3V)/2
                    #                                            = 1 (mod 3^threeadicity)
                    #    343*2^{6j}*X^6   + 74*2^{3j}*X^3   + 3  = 0 (mod ")
                    #
                    solveSextic(343*2^(6*j), 74*2^(3*j), 3, pow3)
            ),
            (17, 1, #    343*x^6          + 17*x^3          + 1  = q + 1 + (T-3V)/2
                    #                                            = 1 (mod 3^threeadicity)
                    #    343*2^{6j}*X^6   + 17*2^{3j}*X^3   + 1  = 1 (mod ")
                    #   (343*2^{6j}*X^3   + 17*2^{3j})*X^3       = 0 (mod ")
                    #
                    # X = 0 or X^3 = -17 / 343*2^{3j} (mod 3^threeadicity)
                    [0] + [int(X_) for X_ in (-17 / (343*Mod(2^(3*j), pow3))).nth_root(3, all=True)]
            ),
            #(0, 0, # 343*x^6 = q + 1 - (T+3V)/2 is never prime
            #       []
            #),
            (57, 3, #    343*x^6          + 57*x^3          + 3  = q + 1 - (T-3V)/2
                    #                                            = 1 (mod 3^threeadicity)
                    # 2*(343*2^{6j-1}*X^6 + 57*2^{3j-1}*X^3 + 1) = 0 (mod ")
                    #    343*2^{6j-1}*X^6 + 57*2^{3j-1}*X^3 + 1  = 0 (mod ")
                    #
                    solveSextic(343*2^(6*j-1), 57*2^(3*j-1), 1, pow3)
            ),
        ]

        def filter_solution(u, v, Xoffsets):
            def check_solution(X):
                x = Mod(X << j, pow3)
                r = int(343*x^6 + u*x^3 + v)
                q = int(343*x^6 + 37*x^3 + 1)
                z = Mod(7*(X << j), pow3*7)
                p21 = int(z^8 + 5*z^7 + 7*z^6 + 37*z^5 + 188*z^4 + 259*z^3 + 343*z^2 + 1763*z + 2401) % 21
                assert r % pow3 == 1, r
                return gcd(r, pow3) == 1 and gcd(q, pow3) == 1 and p21 == 0

            return (u, v, filter(check_solution, Xoffsets))

        solutions = list(starmap(filter_solution, solutions))

        print("Xbase = %s, 2-adicity = %d, 3-adicity = %d" % (format_int(Xbase, 2), twoadicity, threeadicity))
        print("solutions = %r" % (solutions,))

        self.params = (Xbase, Xend, solutions, j, pow2, pow3, L)

    def run(self, wid, processes):
        (Xbase, Xend, solutions, j, pow2, pow3, L) = self.params

        # Align Xbase for this worker.
        chunksize = pow3*processes
        Xbase = ((Xbase+chunksize-1) // chunksize)*chunksize + wid*pow3

        for Xchunk in range(Xbase, Xend, chunksize):
            for (u, v, Xoffsets) in solutions:
                rdesc = "343*x^6 + %d*x^3 + %d" % (u, v)
                for Xoffset in Xoffsets:
                    x = (Xchunk + Xoffset) << j

                    q = 343*x^6 + 37*x^3 + 1
                    if q < 2^L or q % pow2 != 1:
                        print("q = %s" % (format_int(q, 2),))
                        continue

                    r = 343*x^6 + u*x^3 + v
                    if r % pow3 != 1:
                        print("r = %s" % (format_int(r, 3),))
                        continue

                    z = 7*x
                    p21 = z^8 + 5*z^7 + 7*z^6 + 37*z^5 + 188*z^4 + 259*z^3 + 343*z^2 + 1763*z + 2401
                    if p21 % 21 != 0:
                        print("p = %s" % (format_int(p),))
                        continue
                    p = p21//21

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

        cofactor3 = 2401*x^2 + 1715*x + 343
        assert(cofactor3 % 3 == 0)
        cofactor = cofactor3//3
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
        if TWIST_SECURITY > 0 and twsecq < TWIST_SECURITY: continue
        (twsecr, twembedr) = twist_security(r, q)
        if TWIST_SECURITY > 0 and twsecr < TWIST_SECURITY: continue

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
        twembeddivq = None if TWIST_SECURITY == 0 else (2*q + 1 - r)/twembedq
        twembeddivr = None if TWIST_SECURITY == 0 else (2*r + 1 - q)/twembedr

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
    if TWIST_SECURITY == 0:
        return (None, None)
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

        if TWIST_SECURITY > 0:
            output += "Eq twist Pollard security = %.1f, embedding degree = (2q + 1 - r)/%d\n" % (twsecq, twembeddivq)
            output += "Er twist Pollard security = %.1f, embedding degree = (2r + 1 - q)/%d\n" % (twsecr, twembeddivr)

        print(output)  # one syscall to minimize tearing

main()
