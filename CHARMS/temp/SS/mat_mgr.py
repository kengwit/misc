from famuls import CHECK

class MAT_MGR:

    _instance_count = 0
    
    def __init__(self,mgr):
        self._mgr = mgr
        CHECK (self._mgr != None ) # , EXCEPTION_NULL_PTR,;);
        self._instance_count += 1
        CHECK (self._instance_count == 1) #, EXCEPTION_ILLEGAL_USE,;);
        self._mm = {} # material map
  
