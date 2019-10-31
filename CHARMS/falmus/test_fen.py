from fen import FEN

for i in range(1,10):
    
    n = FEN(id=1,ref_loc=0,gmesh=0)
    print(n.uniqobjid())
    assert(n.uniqobjid()==(i-1))
    

