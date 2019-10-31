from famuls import CHECK
from gmesh import GMESH

class MESH_MGR:

    _instance_count = 0
    
    def __init__(self,mgr):
        self._mgr = mgr
        CHECK (self._mgr != None ) # , EXCEPTION_NULL_PTR,;);
        self._instance_count += 1
        CHECK (self._instance_count == 1) #, EXCEPTION_ILLEGAL_USE,;);
        self._mesh_map = {}
        
        
    def gmesh(self,name):
        if ( name in self._mesh_map ):
            return self._mesh_map[name]
        else:            
            gmesh = GMESH(name,self._mgr.db())
            CHECK(gmesh != None) # , EXCEPTION_NULL_PTR,;);
            self._mesh_map[name] = gmesh
            return gmesh
        

class MAT_MGR:

    _instance_count = 0
    
    def __init__(self,mgr):
        self._mgr = mgr
        CHECK (self._mgr != None ) # , EXCEPTION_NULL_PTR,;);
        self._instance_count += 1
        CHECK (self._instance_count == 1) #, EXCEPTION_ILLEGAL_USE,;);
        self._mm = {} # material map
        
        
class LOAD_MGR:

    _instance_count = 0
    
    def __init__(self,mgr):
        self._mgr = mgr
        CHECK (self._mgr != None ) # , EXCEPTION_NULL_PTR,;);
        self._instance_count += 1
        CHECK (self._instance_count == 1) #, EXCEPTION_ILLEGAL_USE,;);
        
        
        
class MGR:
    
    _instance_count = 0
      
    def __init__(self,db):
        self._db = db
        CHECK (self._db != None ) # , EXCEPTION_NULL_PTR,;);
        self._instance_count += 1
        CHECK (self._instance_count == 1) #, EXCEPTION_ILLEGAL_USE,;);
        
        self._mesh_mgr = MESH_MGR(self)
        self._mat_mgr  = MAT_MGR(self)
        self._load_mgr = LOAD_MGR(self)


    # Get the database handle.
    def db(self):
        return self._db
    
    def mesh_mgr(self):
        return self._mesh_mgr
    
    
