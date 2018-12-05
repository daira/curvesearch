import sys
from math import log
from multiprocessing import Pool
from random import randint
from time import sleep
from traceback import print_exc
import os


FIELD_SIZE = 762
#FIELD_SIZE = 2048
TWO_ADICITY = 24
#TWO_ADICITY = 1
CM_THRESHOLD = 10^16
#DETERMINISTIC = False
DETERMINISTIC = True

# tuning
PROCESSES = None  # auto-detect
TIMEOUT = 20

INNER_LOOP = 100
PRIME_TRIALS = 10000


# Consider a "virtual factor wheel" of circumference
# 2^TWO_ADICITY * 3^2 * 5 * 7 * 11 * 13 * 17.
# (See <https://en.wikipedia.org/wiki/Wheel_factorization> for a
# description of factor wheels in general.)
#
# Since we only choose primes of the form k*2^TWO_ADICITY + 1,
# we can compress this wheel to size COMPRESSED_CIRCUMFERENCE:

COMPRESSED_CIRCUMFERENCE = 3^2 * 5 * 7 * 11 * 13 * 17
VIRTUAL_CIRCUMFERENCE = 2^TWO_ADICITY * COMPRESSED_CIRCUMFERENCE
#
# The entry at each index i gives the skip distance from
# (i*2^TWO_ADICITY + 1) to the next integer coprime to
# VIRTUAL_CIRCUMFERENCE. So when searching for primes q2, we can
# efficiently skip between numbers coprime to VIRTUAL_CIRCUMFERENCE,
# which are more likely to be prime.
#
# (Traditionally, a factor wheel would represent only the
# skip distances between consecutive coprime entries, but we also
# need to be able to get from an arbitrary base modulo
# VIRTUAL_CIRCUMFERENCE to the first entry after that base.
# Only representing the coprime entries would have slightly
# better access locality, but that is more complex and actually
# doesn't make much difference for the COMPRESSED_CIRCUMFERENCE
# used here.)
#
# But wait, there's more. For MNT6-BN12 2-cycles, the difference q2 - q1
# is given by 6 * x2^2. Suppose we choose x2 to be a multiple of
#X2_MULTIPLE = 2^(TWO_ADICITY//2) * 3 * 5 * 7 * 11 * 13 * 17

if DETERMINISTIC:
    X2_MULTIPLE = 2^180
else:
    X2_MULTIPLE = 2^166 * 3 * 5 * 7 * 11 * 13 * 17

#
# Then, 6 * x2^2 will be a multiple of
# 2^TWO_ADICITY * 3^3 * 5^2 * 7^2 * 11^2 * 13^2 * 17^2.
#
# Since q2 is prime, it is coprime to this multiple, and so therefore
# will be q1. This makes q1 also more likely to be prime.

TWOADIC = 2^TWO_ADICITY
WHEEL = [0]*COMPRESSED_CIRCUMFERENCE

# We construct the wheel backwards, so that we don't need to look ahead.
skip_distance = 2*TWOADIC
for i in xrange(COMPRESSED_CIRCUMFERENCE-1, -1, -1):
    WHEEL[i] = skip_distance

    if gcd(i*TWOADIC + 1, VIRTUAL_CIRCUMFERENCE) == 1:
        skip_distance = 0

    skip_distance += TWOADIC



def find_mnt6bn(field_size):
    if field_size < 762:
        print("*** Field size %d bits is less than the recommended minimum "
              "of 762 bits. You have been warned. ***\n" % (field_size,))

    processes = PROCESSES or sage.parallel.ncpus.ncpus()
    print("Using %d processes." % (processes,))
    pool = Pool(processes=processes)

    try:
        for w in xrange(processes):
            pool.apply_async(worker, (field_size, w, processes))

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


def real_worker(field_size, w, processes):
    if DETERMINISTIC:
        t2_minus_1 = int(sqrt((2^field_size) // 4))
        x2 = ((int(sqrt(t2_minus_1//6)) + X2_MULTIPLE - 1)//X2_MULTIPLE + w - processes)*X2_MULTIPLE

    while True:
        sys.stdout.write("%d" % (w,))
        sys.stdout.flush()
        if DETERMINISTIC:
            x2 += processes*X2_MULTIPLE
        else:
            t2 = 2^(field_size//2 + 1)
            t2 = randint(t2, (7*t2)//5)
            x2 = int(sqrt(t2//6))
            x2 = ((x2 + X2_MULTIPLE - 1)//X2_MULTIPLE)*X2_MULTIPLE

#        t2 = 6*x2^2 + 1
        x1 = 3*x2^2
        t1 = -2*x1 + 1

        base = t1^2 // 4
        for i in xrange(INNER_LOOP):
            q1 = next_pseudoprime_onemodn(base)
            if q1 is None:
                break

            base = q1  # for next iteration
            sys.stdout.write(".")
            sys.stdout.flush()

            tt = t1^2
            assert(tt < 4*q1)

            #q1 = q2 - t2 + 1
            q2 = q1 - t1 + 1
            if is_pseudoprime(q2):
                D = squarefree_part(4*q1 - tt)
                if D < CM_THRESHOLD and is_prime(q1) and is_prime(q2):
                    print("\nq1 = 0x%x = 2^%f;\nq2 = 0x%x = 2^%f;\nx2 = 0x%x; D = %s = 10^%f" %
                          (q1, log(q1)/log(2), q2, log(q2)/log(2), x2, D, log(D)/log(10)))
                    show_factors("q1-1", q1-1)
                    show_factors("q2-1", q2-1)
                else:
                    sys.stdout.write("#%d" % (log(D)/log(10),))
                    sys.stdout.flush()


@fork(timeout=TIMEOUT)
def show_factors(s, n):
    print("\n0x%x (%s) = %s" % (n, s, factor(n)))


def next_pseudoprime_onemodn(base):
    """
    Return the next pseudoprime *greater* than base that is 1 modulo TWOADIC.
    """

    x = ((base+TWOADIC-1)//TWOADIC + 1)*TWOADIC + 1

    for trials in xrange(PRIME_TRIALS):
        x += WHEEL[((x % VIRTUAL_CIRCUMFERENCE)-1)//TWOADIC]
        if is_pseudoprime(x):
            return x

    return None


find_mnt6bn(FIELD_SIZE)
