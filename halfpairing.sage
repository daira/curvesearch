#!/usr/bin/env sage
# -*- coding: utf-8 -*-

import sys
from multiprocessing import Pool, cpu_count
from traceback import print_exc
from math import ceil
from itertools import combinations

if sys.version_info[0] == 2:
    range = xrange


# Let Ep/Fp : y^2 = x^3 + bp
# Let Eq/Fq : y^2 = x^3 + bq

# p and q should each be ~ L bits.

DEFAULT_TWOADICITY = 32

ALLOW_NEGATIVE_U = True

#COEFFICIENT_RANGE = (5,)
# TODO: could optimize by observing that b can't be square in the field.
# (An elliptic curve of the form y^2 = x^3 + b over F_p where b is square in F_p and
# p is a large prime, always has order a multiple of 3. See Theorems 5.3 and 5.1c
# in 'Elliptic Curves' by Anthony W. Knapp.)
COEFFICIENT_RANGE = range(1, 100)

GCD_PRIMES = (5, 7, 11, 13, 17)

# Set to a prime, or 0 to disable searching for isogenies.
#ISOGENY_DEGREE_MAX = 3
ISOGENY_DEGREE_MAX = 13

DEFAULT_TWIST_SECURITY = 120
REQUIRE_PRIMITIVE = False #True


def low_hamming_order(L, twoadicity, wid, processes):
    # p-1 is a multiple of 6u, so will have exactly one more trailing zero than u when u is odd.
    trailing_zeros = twoadicity-1

    # If u = 2^l, then p ~ 36 * 2^{4l}.
    assert L >= 22
    l = (L - 18)//4
    base = ((8, 10, 11, 14)[(L - 18) % 4] << l) + (1 << trailing_zeros)
    for w in range(wid, l-trailing_zeros+1, processes):
        for c in combinations(range(trailing_zeros+1, l), w):
            u = base + sum([1 << i for i in c])
            yield u
            if ALLOW_NEGATIVE_U:
                yield -u

def find_bn_primes(L, twoadicity, wid, processes):
    for u in low_hamming_order(L, twoadicity, wid, processes):
        p = 36*u^4 + 36*u^3 + 24*u^2 + 6*u + 1
        assert p % (1<<twoadicity) == 1

        if is_pseudoprime(p):
            sys.stdout.write('.')
            sys.stdout.flush()

            q = p - 6*u^2
            if is_prime(q) and is_prime(p):
                assert q % (1<<twoadicity) == 1
                yield (p, q, u)


def find_nice_curves(L, twoadicity, requireisos, sortpq, sortqp, twistsec, wid, processes):
    for (p, q, u) in find_bn_primes(L, twoadicity, wid, processes):
        if (sortpq and q < p) or (sortqp and p < q):
            (p, q) = (q, p)

        (Ep, bp) = find_curve(p, q)
        if bp is None: continue
        #(Eq, bq) = find_curve(q, p, (bp,))
        (Eq, bq) = find_curve(q, p)
        if bq is None: continue

        sys.stdout.write('*')
        sys.stdout.flush()

        primp = (Mod(bp, p).multiplicative_order() == p-1)
        if REQUIRE_PRIMITIVE and not primp: continue
        primq = (Mod(bq, q).multiplicative_order() == q-1)
        if REQUIRE_PRIMITIVE and not primq: continue

        (twsecp, twembedp) = twist_security(p, q)
        if twsecp < twistsec: continue
        (twsecq, twembedq) = twist_security(q, p)
        if twsecq < twistsec: continue

        (secp, embedp) = curve_security(p, q)
        (secq, embedq) = curve_security(q, p)

        zetap = GF(p).zeta(3)
        zetap = min(zetap, zetap^2)
        assert zetap^3 == Mod(1, p)

        zetaq = GF(q).zeta(3)
        P = Ep.gens()[0]
        zP = endo(Ep, zetap, P)
        if zP != int(zetaq)*P:
            zetaq = zetaq^2
            assert(zP == int(zetaq)*P)
        assert zetaq^3 == Mod(1, q)

        Q = Eq.gens()[0]
        assert endo(Eq, zetaq, Q) == int(zetap)*Q

        iso_Ep = find_iso(Ep)
        iso_Eq = find_iso(Eq)
        if requireisos and (iso_Ep is None or iso_Eq is None):
            continue

        yield (p, q, u, bp, bq, zetap, zetaq, primp, primq, secp, secq, twsecp, twsecq,
               embedp, embedq, twembedp, twembedq, iso_Ep, iso_Eq)

def endo(E, zeta, P):
   (xp, yp) = P.xy()
   return E(zeta*xp, yp)

def find_curve(p, q, b_range=None):
    for b in (b_range or COEFFICIENT_RANGE):
        E = EllipticCurve(GF(p), [0, b])
        if E.count_points() == q:
            return (E, b)
    return (None, None)

def find_gcd_primes(p):
    return (r for r in GCD_PRIMES if gcd(p-1, r) == 1)

pi_12 = (pi/12).numerical_approx()

def curve_security(p, q):
    sys.stdout.write('!')
    sys.stdout.flush()
    r = factor(q)[-1][0]
    return (log(pi_12 * r, 4), embedding_degree(p, r))

def twist_security(p, q):
    return curve_security(p, 2*(p+1) - q)

def embedding_degree(p, r):
    sys.stdout.write('#')
    sys.stdout.flush()
    assert gcd(p, r) == 1
    Z_q = Integers(r)
    u = Z_q(p)
    d = r-1
    V = factor(d)
    for (v, k) in V:
        while d % v == 0:
            if u^(d/v) != 1: break
            d /= v

    return d

def find_iso(E):
    # Based on <https://eprint.iacr.org/2019/403.pdf> Appendix A.
    # Look for isogenous curves having j-invariant not in {0, 1728}.
    for degree in primes(ISOGENY_DEGREE_MAX+1):
        sys.stdout.write('~')
        sys.stdout.flush()
        for iso in E.isogenies_prime_degree(degree):
            if iso.codomain().j_invariant() not in (0, 1728):
                return iso.dual()

    return None


def format_weight(x, detail=True, weight=False):
    X = format(abs(x), 'b')
    if detail:
        assert X.endswith('1')
        detailstr = " (bitlength %d, weight %d, 2-adicity %d)" % (len(X), sum([int(c) for c in X]),
                                                                  len(X) - len(X[:-1].rstrip('0')))
    elif weight:
        detailstr = " (bitlength %d, weight %d)" % (len(X), sum([int(c) for c in X]))
    else:
        detailstr = " (bitlength %d)" % (len(X),)

    return "%s0b%s%s" % ("-" if x < 0 else "", X, detailstr)

def format_embedding_degree(ndesc, n, embed):
    if embed < 100: return "%d" % (embed,)
    return "(%s)/%d" % (ndesc, n/embed)


def main():
    args = sys.argv[1:]
    processes = 1 if "--sequential" in args else cpu_count()
    if processes >= 6:
        processes -= 2
    requireisos = "--requireisos" in args
    sortpq = "--sortpq" in args
    sortqp = "--sortqp" in args
    twistsec = 0 if "--ignoretwist" in args else DEFAULT_TWIST_SECURITY
    args = [arg for arg in args if not arg.startswith("--")]

    if len(args) < 1:
        print("""
Usage: sage halfpairing.sage [--sequential] [--requireisos] [--sortpq] [--sortqp] [--ignoretwist]
                             <min-bitlength> [<min-2adicity>]

Arguments:
  --sequential     Use only one thread, avoiding non-determinism in the output order.
  --requireisos    Require isogenies useful for a "simplified SWU" hash to curve.
  --sortpq         Sort p smaller than q.
  --sortqp         Sort q smaller than p.
  --ignoretwist    Ignore twist security.
  <min-bitlength>  Both primes should have this minimum bit length.
  <min-2adicity>   Both primes should have this minimum 2-adicity.
""")
        return

    L          = int(args[0])
    twoadicity = int(args[1]) if len(args) > 1 else DEFAULT_TWOADICITY

    print("Using %d processes." % (processes,))
    pool = Pool(processes=processes)

    try:
        for wid in range(processes):
            pool.apply_async(worker, (L, twoadicity, requireisos, sortpq, sortqp, twistsec, wid, processes))

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
    for (p, q, u, bp, bq, zetap, zetaq, primp, primq, secp, secq, twsecp, twsecq,
         embedp, embedq, twembedp, twembedq, iso_Ep, iso_Eq) in find_nice_curves(*args):
        output  = "\n"
        output += "p   = %s\n" % format_weight(p)
        output += "q   = %s\n" % format_weight(q)
        output += "u   = %s\n" % format_weight(u, detail=False, weight=True)
        output += "ζ_p = %s (mod p)\n" % format_weight(int(zetap), detail=False)
        output += "ζ_q = %s (mod q)\n" % format_weight(int(zetaq), detail=False)

        output += "Ep/Fp : y^2 = x^3 + %d\n" % (bp,)
        output += "Eq/Fq : y^2 = x^3 + %d\n" % (bq,)

        output += "gcd(p-1, α) = 1 for α ∊ {%s}\n" % (", ".join(map(str, find_gcd_primes(p))),)
        output += "gcd(q-1, α) = 1 for α ∊ {%s}\n" % (", ".join(map(str, find_gcd_primes(q))),)

        output += "%d is %sprimitive in Fp\n" % (bp, "" if primp else "non")
        output += "%d is %sprimitive in Fq\n" % (bq, "" if primq else "non")

        output += "Ep rho security = %.1f, embedding degree = %s\n" % (secp, format_embedding_degree("q-1", q-1, embedp))
        output += "Eq rho security = %.1f, embedding degree = %s\n" % (secq, format_embedding_degree("p-1", p-1, embedq))

        output += "Ep twist rho security = %.1f, embedding degree = %s\n" % (twsecp, format_embedding_degree("2p + 1 - q", 2*p + 1 - q, twembedp))
        output += "Eq twist rho security = %.1f, embedding degree = %s\n" % (twsecq, format_embedding_degree("2q + 1 - p", 2*q + 1 - p, twembedq))

        if iso_Ep is not None:
            output += "iso_Ep = %r\n" % (iso_Ep,)
            output += "iso_Ep maps = %r\n" % (iso_Ep.rational_maps(),)
        elif ISOGENY_DEGREE_MAX > 0:
            output += "No Ep isogenies for simplified SWU with degree ≤ %d\n" % (ISOGENY_DEGREE_MAX,)

        if iso_Eq is not None:
            output += "iso_Eq = %r\n" % (iso_Eq,)
            output += "iso_Eq maps = %r\n" % (iso_Eq.rational_maps(),)
        elif ISOGENY_DEGREE_MAX > 0:
            output += "No Eq isogenies for simplified SWU with degree ≤ %d\n" % (ISOGENY_DEGREE_MAX,)

        print(output)  # one syscall to minimize tearing

main()
