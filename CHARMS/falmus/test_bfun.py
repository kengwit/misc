from bfun import *
from gmesh import GMESH, GSUBMESH
from bfun import BFUN_SET

from algo import ALGO_ELASTICITY
from db import DB
from mgr import MGR



jobname = 'Mandible'

db = DB()
db.load_file(jobname + '.yaml')

mgr = MGR(db)

a = ALGO_ELASTICITY(jobname, mgr);
a.setup()

gmesh = a.gmesh()
gsubmesh = a._gsubmeshes[0]

bfun_set = BFUN_SET.make(BFUN_SET.REFINEMENT_STRATEGY.QUASI_HIERARCHICAL)