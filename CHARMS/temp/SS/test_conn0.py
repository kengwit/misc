from conn import *


typedef_CL2 = CONN_LINE_2

c1 = typedef_CL2()
#cloned_c1 = c1.clone()
#
#assert(c1.MANIFOLD_DIM == cloned_c1.MANIFOLD_DIM)
#assert(c1.NFENS == cloned_c1.NFENS)
#assert(c1.NREFFENS == cloned_c1.NREFFENS)
#
#
#CL2a = CONN_LINE_2()
#CL2b = CONN_LINE_2()
#CL2a._fens=[1,2]
#CL2b._fens=[1,2]
#
#cloned_CL2a = CL2a.clone()
#print(cloned_CL2a._fens)
#
#ta = CL2a.get_ref_fen(CL2b,ref_fens=[3],indx=0)
#cloned_ta = cloned_CL2a.get_ref_fen(CL2b,ref_fens=[3],indx=0)