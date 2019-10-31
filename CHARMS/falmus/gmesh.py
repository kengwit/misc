from famuls import CHECK
from fen import FEN
from point import POINT
from gcell import GCELL_GROUP
  


class GSUBMESH:
    
    class gcell_group_enumerator_t():        
        
        def __init__(self,gcell_groups):
            self._gcell_groups = gcell_groups
            
        def __iter__(self):
            for gcell_group in self._gcell_groups:
                yield gcell_group
                
    def gcell_group_enumerator(self):
        e = self.gcell_group_enumerator_t(self._gcell_groups)
        return e
    
    def __init__(self,name,db,gmesh):
        self._name = name
        self._gcell_groups = []
        self._gmesh = gmesh
        
        # GCELL_GROUP
        path = "algorithms/gsubmeshes/" + self._name + "/gcell_groups"
        gcell_group_list = db.DB_GET_STRING_LIST(path)
        CHECK ( len( gcell_group_list ) != 0 ) # EXCEPTION_BAD_VALUE,;);
        for gg_name in gcell_group_list:
            gcell_group = GCELL_GROUP(gg_name,db,self)
            self._gcell_groups.append(gcell_group)

    def gmesh(self):
        return self._gmesh        
    
    def name(self):
        return self._name
        
class GMESH:
    
    class gsubmesh_enumerator_t():        
        
        def __init__(self,gsubmeshes):
            self._gsubmeshes = gsubmeshes
            
        def __iter__(self):
            for gsubmesh in self._gsubmeshes:
                yield gsubmesh
    
    class fens_enumerator_t():        
        
        def __init__(self,fens):
            self._fens = fens
            
        def __iter__(self):
            for fen in self._fens:
                yield fen
    
    # Get an enumerator of the component submeshes.
    def gsubmesh_enumerator(self):
        e = self.gsubmesh_enumerator_t(self._gsubmeshes)
        return e
    
    # Get an enumerator of the finite element nodes.
    def fens_enumerator(self):
        e = self.fens_enumerator_t(self._fens)
        return e
    
    
    def __init__(self,name,db):
        self._name = name
        self._fens       = []
        self._gsubmeshes = []
        self._fen_map    = {}
        self._max_fen_id = 0

        # FEN
        fens_param = db.DB_GET_STRING ('algorithms/gmeshes/' + self._name + '/fen_file')
        CHECK (fens_param != None) #, EXCEPTION_BAD_VALUE,;);
        CHECK (self.read_fens(fens_param) ) #, EXCEPTION_BAD_VALUE,;);
        
        # GSUBMESH        
        path = "algorithms/gmeshes/" + self._name + "/gsubmeshes"
        gsubmeshes_list = db.DB_GET_STRING_LIST(path)
        
        CHECK (len(gsubmeshes_list) != 0 ) #, EXCEPTION_BAD_VALUE,;);
        for gsm in gsubmeshes_list:
            gsubmesh = GSUBMESH(gsm,db,self)
            self._gsubmeshes.append(gsubmesh)
        
        
        
    def name(self):
        return self._name
    
 
    # Get the total of nodes associated with mesh.
    def fen_total(self):
        return len(self._fens)
    
        
    # Find a finite element node by id.  
    def find_fen(self,id):
        if (id in self._fen_map):
            return self._fen_map[id]
        else:
            return None
    
    def find_gsubmesh(self,name):
        e = self.gsubmesh_enumerator()      
        #for gsm in self._gsubmeshes:
        for gsm in e:
            if ( gsm.name() == name ):
                return gsm
        
        return None # no such submesh

    # Add a finite element node.
    def add_fen(self,fen):
        
        # Generate an internal id for the node
        # Add to a list
        self._fens.append(fen)
        
        # Largest external node id (to be used during refinement)
        self._max_fen_id = max(self._max_fen_id, fen.id())
        
    
    def max_fen_id (self):
        return self._max_fen_id
        
    def read_fens(self,file):
        
        read_flag = False
        content=[]
        with open(file,'r') as f:
            for line in f:
                content.append(line)
            
        for line in content:
            line = line.strip().split()
            id = int(line[0])
            p = POINT()
            p[0] = float(line[1])
            p[1] = float(line[2])
            p[2] = float(line[3])
            fen = FEN(id,p,self)
            self._fen_map[fen.id()] = fen
            
        read_flag = True
        
        return read_flag
    