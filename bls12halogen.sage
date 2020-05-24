# -*- coding: utf-8 -*-
import sys
from multiprocessing import Pool, cpu_count
from traceback import print_exc
from math import ceil
from itertools import combinations

if sys.version_info[0] == 2: range = xrange

# y^2 = x^3 + 1 never has amicable pairs. <https://arxiv.org/pdf/0912.1831.pdf>
COEFFICIENT_RANGE_Q = range(2, 100)
COEFFICIENT_RANGE_R = range(2, 100)

DEFAULT_TWOADICITY = 32


class BruteForce:
    def __init__(self, L, twoadicity):
        j = (twoadicity+1)//2
        threeadicity = ceil(twoadicity*log(2, 3))

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
        # This may or may not have a solution for X. If it does, X is heuristically
        # about the same length as 3^threeadicity, and so we have approximately
        # A ~ L/4 - 2j bits.

        Alen = (L+3)//4 - 2*j
        Abase = 1 << Alen
        Aend = Abase * 4
        self.pow2 = 2^twoadicity
        self.pow3 = 3^threeadicity
        self.params = (Abase, Aend, j, L)

    def run(self, wid, processes):
        (Abase, Aend, j, L) = self.params

        # Align Abase for this worker.
        Abase = ((Abase+processes-1) // processes)*processes + wid

        for A in range(Abase, Aend, processes):
            X4 = -1 / Mod(2^(4*j - 1) * A^4, self.pow3)
            if not sqrt(X4).is_square():
                continue
            X = int(sqrt(sqrt(X4)))
            x = (A * X) << j

            # p is less likely to be prime than r, so check p first.
            #print("A = %s, X = %s, x = %s = %s" % (format_int(A), format_int(X, 2), format_int(x, 2), format_int(x, 3)))
            if x % 3 == 1:
                q = x^4 - x^2 + 1
                #print("q = %s" % (format_int(q, 2),))
                if q < 2^L:
                    continue
                assert(q % self.pow2 == 1)
                assert((q*(x-1)^2) % 3 == 0)
                p = x + (q*(x-1)^2)//3
                #print("p = %s" % (format_int(p, 2),))
                if is_pseudoprime(p):
                    sys.stderr.write('.')
                    sys.stderr.flush()
                    if is_pseudoprime(q):
                        yield (p, q, x)

        sys.stderr.write('<')
        sys.stderr.flush()

def find_cycle(q, pow2, pow3):
    set_r = set()
    for bq in COEFFICIENT_RANGE_Q:
        Eq = EllipticCurve(GF(q), [0, bq])
        r = Eq.count_points()
        assert((r - (r+1))^2 < 4*q)
        t = q + 1 - r
        #print("\nq = %s, r = %s, trace = %s, r = %d (mod 126)" % (format_int(q, 2), format_int(r, 3), format_int(t), int(r)%126))
        #if r % pow3 == 1 and:
        if is_pseudoprime(r):
            #print("\n%4d: %s (%d, %d, %2d (mod 126))" % (bq, format_int(r), int(q).bit_length(), int(r).bit_length(), int(r)%126))
            sys.stdout.write('#')
            sys.stdout.flush()
            set_qdash = set()
            for br in COEFFICIENT_RANGE_R:
                sys.stdout.write('!')
                sys.stdout.flush()
                Er = EllipticCurve(GF(r), [0, br])
                qdash = Er.count_points()
                if qdash == q:
                    return (Eq, bq, Er, br, r)

                set_qdash.add(qdash)
                if len(set_qdash) == 6: break

        set_r.add(r)
        if len(set_r) == 6: break

    return (None, None, None, None, None)

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

def format_int(n, b=10):
    if n is None: return 'None'
    n = int(n)
    if n == 0: return '0'
    neg = " " if n > 0 else "-"
    n = abs(n)
    nums = []
    while n > 0:
        n, r = divmod(n, b)
        nums.append(str(r))
    return "%s%s_%d" % (neg, ''.join(reversed(nums)), b)

def real_worker(*args):
    (strategy, wid, processes) = args
    for (p, q, x) in strategy.run(wid, processes):
        (Eq, bq, Er, br, r) = find_cycle(q, strategy.pow2, strategy.pow3)
        if Er is None: continue
        output  = "\n"
        output += "p  = %s\n" % format_int(p)
        output += "q  = %s\n" % format_int(q, 2)
        output += "r  = %s\n" % format_int(r, 3)
        output += "x  = %s\n" % format_int(x, 2)
        output += "Eq = %s\n" % Eq
        output += "Er = %s\n" % Er

        print(output)  # one syscall to minimize tearing

main()
