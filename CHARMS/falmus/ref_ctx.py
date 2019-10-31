from ref_fen_src import REF_FEN_SRC
from fen_maker import FEN_MAKER

class REF_CTX:
    '''
     A refinement context.
     A refinement context holds the incremental changes to the basis
     function set resulting from refinement or derefinement.

     A refinement context holds a set of active
     nodes and nodes to be deactivated.  It takes 
     a field (or perhaps several fields) and refines the geometric
     cells to build up a geometrical structure for basis function refinement.
     The nodes of the basis functions to deactivate are also listed
     in the refinement context as inactive nodes.  
    '''
    
    def __init__(self,gmesh):
        self._gmesh = gmesh
        self._ref_fen_src = REF_FEN_SRC()
        self._activated_fens = set([])
        self._deactivated_fens = set([])
        self._true_hierarchical = True
        
    def true_hierarchical(self):
        return self._true_hierarchical
    
    
    def set_true_hierarchical(self,true_hierarchical):
        self._true_hierarchical = true_hierarchical
        
    
    
    # Get a refinement node --- see REF_FEN_SRC.
    def get_ref_fen (self, conn, indx, ref_locs):        
        fen_maker = FEN_MAKER(self._gmesh, ref_locs)
        
        raise NotImplementedError
        ref_fen_maker = [fen_maker] # wrap fen_maker with a list so that it becomes a reference
        self._ref_fen_src.get_ref_fen(conn, indx, ref_fen_maker)
        fen_maker = ref_fen_maker[0]