import sys
from itertools import chain

JUST_PRINT = True

x = QQ['x'].0

def run():
    ss = []
    process_candidates = add_candidates if JUST_PRINT else check_for_2cycles

    # u(x) must be always-even, so that n and p can both be prime.

    for b in xrange(0, 20, 2):
        for c in xrange(0, b, 2):
            process_candidates(ss, 0, b, c)

    for a in chain(xrange(-20, 0), xrange(1, 20)):
        # For any integer z we can substitute x -> x+z giving a*(x+z)^2 = a*x^2 + 2*a*z*x + a*z^2.
        # Also, negative b is equivalent to positive b under the mapping x -> -x.
        # Therefore, we need not consider b outside [0, 2*a).
        for b in xrange(0, 2*a):
            # If a and b are both even or both odd, then a*x^2 + b*x will be always-even.
            if (a+b) % 2 == 0:
                crange = xrange(-20, 20, 2)
            else:
                crange = xrange(-20, 20)

            for c in crange:
                process_candidates(ss, a, b, c)

    if JUST_PRINT:
        sys.stderr.write('\n')
        ss.sort(key=lambda (k, t, n, p, DV2, D): (k, D.degree(), t, D))
        print "\n".join(["k=%2s; t=%20s; n=%s; p=%s; DV^2=%s; D=%s" % s for s in ss])


def add_candidates(ss, a, b, c):
    u = a*x^2 + b*x + c
    Phi3  = sum([u^i for i in xrange(0, 3)])
    Phi4  = u^2 + 1
    Phi5  = sum([u^i for i in xrange(0, 5)])
    Phi6  = u^2 - u + 1
    Phi7  = sum([u^i for i in xrange(0, 7)])
    Phi8  = u^4 + 1
    Phi9  = u^6 + u^3 + 1
    Phi10 = u^4 - u^3 + u^2 - u + 1
    Phi11 = sum([u^i for i in xrange(0, 11)])
    Phi12 = u^4 - u^2 + 1
    Phi13 = sum([u^i for i in xrange(0, 13)])
    Phi14 = u^6 - u^5 + u^4 - u^3 + u^2 - u + 1
    Phi15 = u^8 - u^7 + u^5 - u^4 + u^3 - u + 1
    Phi16 = u^8 + 1
    Phi17 = sum([u^i for i in xrange(0, 17)])
    Phi18 = u^6 - u^3 + 1
    Phi19 = sum([u^i for i in xrange(0, 19)])
    Phi20 = u^8 - u^6 + u^4 - u^2 + 1
    Phi21 = u^12 - u^11 + u^9 - u^8 + u^6 - u^4 + u^3 - u + 1
    Phi22 = u^10 - u^9 + u^8 - u^7 + u^6 - u^5 + u^4 - u^3 + u^2 - u + 1
    Phi23 = sum([u^i for i in xrange(0, 23)])
    Phi24 = u^8 - u^4 + 1
    Phi25 = u^20 + u^15 + u^10 + u^5 + 1
    Phi26 = u^12 - u^11 + u^10 - u^9 + u^8 - u^7 + u^6 - u^5 + u^4 - u^3 + u^2 - u + 1
    Phi27 = u^18 + u^9 + 1
    Phi28 = u^12 - u^10 + u^8 - u^6 + u^4 - u^2 + 1
    Phi29 = sum([u^i for i in xrange(0, 29)])
    Phi30 = u^8 - u^4 + 1

    cs = []
    # For cases where the Phi polynomial has no odd powers of u,
    # we needn't consider negative a.
    cs.append(( 3, Phi3))
    if a >= 0: cs.append(( 4, Phi4))
    cs.append(( 5, Phi5))
    cs.append(( 6, Phi6))
    cs.append(( 7, Phi7))
    if a >= 0: cs.append(( 8, Phi8))
    cs.append(( 9, Phi9))
    cs.append((10, Phi10))
    cs.append((11, Phi11))
    if a >= 0: cs.append((12, Phi12))
    cs.append((13, Phi13))
    cs.append((14, Phi14))
    cs.append((15, Phi15))
    if a >= 0: cs.append((16, Phi16))
    cs.append((17, Phi17))
    cs.append((18, Phi18))
    cs.append((19, Phi19))
    if a >= 0: cs.append((20, Phi20))
    cs.append((21, Phi21))
    cs.append((22, Phi22))
    cs.append((23, Phi23))
    if a >= 0: cs.append((24, Phi24))
    cs.append((25, Phi25))
    cs.append((26, Phi26))
    cs.append((27, Phi27))
    if a >= 0: cs.append((28, Phi28))
    cs.append((29, Phi29))
    cs.append((30, Phi30))

    for (k, Phi) in cs:
        factors = integerize(factor(Phi))
        t = u + 1
        for (n, power) in factors:
            p = n + u
            factors_p = integerize(factor(p))
            if factors_p.unit() == 1 and len(factors_p) == 1:
                DV2 = 4*p - t^2
#                print "k=%s; t=%s; factors=%s; n=%s; power=%s; p=%s; DV^2=%s" % (k, t, factors, n, power, p, DV2)
                if (n.degree() == p.degree() and p.coefficients(sparse=False)[0] != 0):
#                    sys.stderr.write('.')
#                    sys.stderr.flush()
                    D = squarefree_part_fixed(DV2)
                    if D.degree() <= 3 and D.coefficients(sparse=False)[D.degree()] > 0:
                        sys.stderr.write(':')
                        sys.stderr.flush()
                        ss.append((k, t, n, p, DV2, D))


def check_for_2cycles(ss, a, b, c):
    ss_1 = []
    ss_2 = []
    add_candidates(ss_1,  a,  b,  c)
    add_candidates(ss_2, -a, -b, -c)
    for (k_1, t_1, n_1, p_1, DV2_1, D_1) in ss_1:
        for (k_2, t_2, n_2, p_2, DV2_2, D_2) in ss_2:
            if n_1 == p_2:
                assert(n_2 == p_1)
                assert(DV2_2 == DV2_1)
                assert(D_2 == D_1)
                print("\nEureka! k_1=%s; k_2=%s; t_1=%s; t_2=%s; n_1=p_2=%s; n_2=p_1=%s; DV2=%s; D=%s" %
                      (k_1, k_2, t_1, t_2, n_1, n_2, DV2_1, D_1))


def integerize(factors):
    assert(isinstance(factors, Factorization))
    u = factors.unit()
    ifactors = []
    for (term, power) in factors:
        denoms = [r.denominator() for r in term.coefficients()]
        scale = lcm(denoms)
        u /= scale^power
        ifactors.append((term*scale, power))

    return Factorization(ifactors, unit=u, simplify=False)


def squarefree_part_fixed(f):
#    print "f=%s" % (f,)
    factors = integerize(factor(f))
#    print "factors=%s" % (factors,)
    u = factors.unit()
    assert(u.denominator() == 1)
    r = squarefree_part(u)
    for (term, power) in factors:
        r *= term^(power % 2)

#    print "r=%s" % (r,)
    return r


run()

