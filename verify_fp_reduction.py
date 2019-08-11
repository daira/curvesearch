#!/bin/env python

from operator import add
from itertools import chain

# Everything here is little-endian.

p = (1<<384)-1

OFFSETS         = [  0, 124, 248]
REDUCED_WIDTHS  = [124, 124, 136]
WIDTHS          = [126, 126, 138]

PRODUCT_OFFSETS = [  0, 248, 372, 496]
PRODUCT_WIDTHS  = [377, 266, 265, 276]

REDUCED_WIDTH_MASKS = [(1<<w)-1 for w in REDUCED_WIDTHS]
PRODUCT_TOTAL_WIDTH = PRODUCT_OFFSETS[-1] + PRODUCT_WIDTHS[-1]

def to_limbs(x):
    return [(x >> OFFSETS[i]) & REDUCED_WIDTH_MASKS[i] for i in range(len(OFFSETS))]

def valid_limbs(limbs):
    return all([limbs[i] < (1 << WIDTHS[i]) for i in range(len(WIDTHS))])

# If all bits were 1, would the reduction fit?
def verify():
    # Compute 2^i mod p, as normalized limbs, for all relevant bit positions.
    indices = sum([range(PRODUCT_OFFSETS[i], PRODUCT_OFFSETS[i]+PRODUCT_WIDTHS[i])
                   for i in range(len(PRODUCT_OFFSETS))], [])

    twoi_modp = [to_limbs(2**i % p) for i in indices]

    # Add them elementwise.
    worst_case = reduce(lambda x, y: map(add, x, y), twoi_modp)
    assert(len(worst_case) == len(OFFSETS))

    print(worst_case)
    print([(1 << WIDTHS[i])-1 for i in range(len(WIDTHS))])
    print(valid_limbs(worst_case))

verify()
