# -*- coding: utf-8 -*-
import sys
from multiprocessing import Pool, cpu_count
from traceback import print_exc
from math import ceil
from itertools import combinations

# Let Ep/Fp : y^2 = x^3 + bp
# Let Eq/Fq : y^2 = x^3 + bq

# p and q should each be ~ L bits.

#COEFFICIENT_RANGE = (5,)
COEFFICIENT_RANGE = xrange(1, 100)

#ACCEPTABLE_PRIMES = (5,)
ACCEPTABLE_PRIMES = Primes()


def find_curve(p):
    for bp in COEFFICIENT_RANGE:
        sys.stdout.write('.')
        sys.stdout.flush()

        for ap in COEFFICIENT_RANGE:
            Ep = EllipticCurve(GF(p), [ap, bp])
            q = Ep.count_points()
            if not is_pseudoprime(q): continue

            for bq in COEFFICIENT_RANGE:
                for aq in COEFFICIENT_RANGE:
                    Eq = EllipticCurve(GF(p), [aq, bq])
                    if Eq.count_points() == p:
                        return (p, ap, bp, q, aq, bq)

    return None

def find_lowest_prime(p):
    for r in ACCEPTABLE_PRIMES:
        if gcd(p-1, r) == 1:
            return r

pi_12 = (pi/12).numerical_approx()

def curve_security(p, q):
    sys.stdout.write('!')
    sys.stdout.flush()
    r = factor(q)[-1][0]
    return (log(pi_12 * r, 4), embedding_degree(p, r))

def embedding_degree(p, r):
    sys.stdout.write('#')
    sys.stdout.flush()
    assert(gcd(p, r) == 1)
    Z_q = Integers(r)
    u = Z_q(p)
    d = r-1
    V = factor(d)
    for (v, k) in V:
        while d % v == 0:
            if u^(d/v) != 1: break
            d /= v

    return d


def format_weight(x, detail=True):
    X = format(abs(x), 'b')
    if detail:
        assert(X.endswith('1'))
        detailstr = " (bitlength %d, weight %d, 2-adicity %d)" % (len(X), sum([int(c) for c in X]),
                                                                  len(X) - len(X[:-1].rstrip('0')))
    else:
        detailstr = " (bitlength %d)" % (len(X),)

    return "%s0b%s%s" % ("-" if x < 0 else "", X, detailstr)


def main():
    p = 52435875175126190479447740508185965837690552500527637822603658699938581184513
    print(find_curve(p))

main()
