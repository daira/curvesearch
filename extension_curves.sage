#!/usr/bin/env sage
# <https://hackmd.io/@daira/HkSXEjD3j>
# Motivated by <https://twitter.com/CPerezz19/status/1619996648838668290>.

from itertools import dropwhile
def red(S): return "\x1b[1m\x1b[31m" + S + "\x1b[22m\x1b[0m"
FACTOR = False

def curve(b, F):
    return EllipticCurve(F, [0, 0, 0, 0, b])

for p in dropwhile(lambda p_: p_ <= 3, Primes()):
    for r in range(2, 20):
        if (p^r) % 3 != 1:
            continue
        F.<u> = GF(p^r)
        B = F.multiplicative_generator()
        candidates = [(i, curve(B^i, F).count_points()) for i in range(1, 7)]
        print("%d^%d = %d, B = %r" % (p, r, p^r, B))
        for (i, s) in candidates:
            print("  %r, %r: %r = %r" % (i, B^i, s, factor(s) if FACTOR else "..."))
            if i in (2, 5): assert s % 2 == 1
            if i in (3, 6): assert s % 2 == 0
            if i in (1, 3, 5): assert s % 3 == 1
            if i in (2, 4, 6): assert s % 3 == 0
            if is_prime(s):
                print(red("  *** FOUND PRIME ***"))
                assert i in (1, 5)

