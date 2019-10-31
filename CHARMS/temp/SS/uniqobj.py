class UNIQOBJ:
    
    _nobj = 0 # static variable
    
    def __init__(self):
        self._internal_id = UNIQOBJ._nobj
        UNIQOBJ._nobj += 1
        
    def id(self):
        return self._internal_id
    
    