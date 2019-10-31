from conn import CONN_MANIFOLD_DIM, CONN_POINT_1, CONN_LINE_2

#typedef_C = CONN(CONN_MANIFOLD_DIM.CONN_0_MANIFOLD
pt   = CONN_POINT_1()
line = CONN_LINE_2()

assert(pt.MANIFOLD_DIM == CONN_MANIFOLD_DIM.CONN_0_MANIFOLD)
assert(line.MANIFOLD_DIM == CONN_MANIFOLD_DIM.CONN_1_MANIFOLD)

pt._fens=[1]
line._fens=[1,2]

cloned_line = line.clone()

pt_ref          = pt.get_ref_fen(pt,ref_fens=[3],indx=0)
line_ref        = line.get_ref_fen(line,ref_fens=[3],indx=0)
cloned_line_ref = cloned_line.get_ref_fen(line,ref_fens=[3],indx=0)
assert(line_ref==cloned_line_ref)