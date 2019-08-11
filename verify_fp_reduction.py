#!/bin/env python

from operator import add
from itertools import chain
from functools import reduce

# Everything here is little-endian.

# This is p for a real 382-bit curve that forms a half-pairing-friendly amicable pair
# (the first one found by `sage halfpairing.sage 382`).
p = 0b1001000000001001000000000011011000000000100100000000000011011000000001011010000000101000100000101011111000011011100111011110110101000100100110101111100100000000100011011011011101100011101001000110110101101000001000001110101000100101000101101000001001011011000110110000000110001100001000011001000011110000000101000100000100100000000001001000000011000000000001100000000000000000000001

OFFSETS         = [  0, 119, 238]
REDUCED_WIDTHS  = [119, 119, 142]
WIDTHS          = [128, 128, 151]

PRODUCT_OFFSETS = [  0, 238, 357, 476]
PRODUCT_WIDTHS  = [377, 304, 303, 302]

REDUCED_WIDTH_MASKS = [(1<<w)-1 for w in REDUCED_WIDTHS]
PRODUCT_TOTAL_WIDTH = PRODUCT_OFFSETS[-1] + PRODUCT_WIDTHS[-1]

def to_limbs(x):
    return [(x >> OFFSETS[i]) & REDUCED_WIDTH_MASKS[i] for i in range(len(OFFSETS))]

def valid_limbs(limbs):
    return all([limbs[i] < (1 << WIDTHS[i]) for i in range(len(WIDTHS))])

# If all bits were 1, would the reduction fit?
def verify():
    # Compute 2^i mod p, as normalized limbs, for all relevant bit positions.
    indices = sum([list(range(PRODUCT_OFFSETS[i], PRODUCT_OFFSETS[i]+PRODUCT_WIDTHS[i]))
                   for i in range(len(PRODUCT_OFFSETS))], [])

    twoi_modp = [to_limbs(2**i % p) for i in indices]

    # Add them elementwise.
    worst_case = list(reduce(lambda x, y: map(add, x, y), twoi_modp))
    assert(len(worst_case) == len(OFFSETS))

    print(worst_case)
    print([(1 << WIDTHS[i])-1 for i in range(len(WIDTHS))])
    print(valid_limbs(worst_case))

verify()
