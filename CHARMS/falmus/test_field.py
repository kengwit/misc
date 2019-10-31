from bfun import BFUN_SET
from field import FIELD_SCALAR, FIELD_VECTOR
from point import FIXED_VECTOR

bfset = BFUN_SET.make(BFUN_SET.REFINEMENT_STRATEGY.QUASI_HIERARCHICAL)

zero = FIXED_VECTOR(3)

fscalar = FIELD_SCALAR("testfield",bfset,zero)
fvector = FIELD_VECTOR("testfield",bfset,zero)

assert(fscalar._bfun_set==bfset)
assert(fvector._bfun_set==bfset)
assert(fscalar.NCOMPONENTS == 1)
assert(fvector.NCOMPONENTS == 3)
