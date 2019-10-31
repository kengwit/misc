class FEN:
    
    # static variable
    _nobj = 0 
	
    def __init__(self, id, ref_loc, gmesh ):
        self._internal_id = FEN._nobj
        FEN._nobj += 1
		
        self._id       = id
        self._ref_loc  = ref_loc
        self._gmesh    = gmesh
        
    def uniqobjid(self):        
        return self._internal_id
			
    # Get the node identifier.
    def id(self):
        return self._id
	