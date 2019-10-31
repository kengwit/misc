from fen import FEN

class FEN_MAKER:
    
    def __init__(self,gmesh,ref_locs):
        self._gmesh = gmesh
        self._ref_locs = ref_locs
    
    def __call__(self, i):
        fen = FEN(self._gmesh.max_fen_id()+1, self._ref_locs[i], self._gmesh)
        return fen
  
