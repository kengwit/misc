from famuls import CHECK
from point import FIXED_VECTOR
from bfun import BFUN_SET
from field import FIELD_VECTOR
from algo import ALGO, ALGO_REFINE

class ALGO_ELASTICITY( ALGO ):
    
    def __init__(self, name, mgr ):
        
        ALGO.__init__(self, name, mgr)
        
        CHECK (self.mgr().db() != None) #, EXCEPTION_NULL_PTR,;);        
        
        self._gmesh       = None        
        top = 'algorithms/elasticity/' + self.name()
        gmesh_name = self.mgr().db().DB_GET_STRING(top + '/gmesh')
        #print(gmesh_name)
        self._gmesh = self.mgr().mesh_mgr().gmesh(gmesh_name)
        
        self._gsubmeshes  = []
        
        self._algo_refine = ALGO_REFINE(name, mgr, self._gmesh)
        self._geometry    = None
        self._u           = None
        self._ebc         = None
        self._algo_errest = None
        
        
    def setup(self):
        CHECK (self.mgr().db() != None) #, EXCEPTION_NULL_PTR,;);
        top = "algorithms/elasticity/" + self.name()
        #print(top)
        # construct the gsubmesh list        
        gsubmesh_names = self.mgr().db().DB_GET_STRING_LIST(top+"/gsubmeshes")
        for gsm in gsubmesh_names:
            gsubmesh = self._gmesh.find_gsubmesh(gsm)
            CHECK(gsubmesh is not None) #, EXCEPTION_NULL_PTR,;);
            self._gsubmeshes.append(gsubmesh);
    
        # create the geometry field
        zero = FIXED_VECTOR(3)
        geom_bfun_set = BFUN_SET.make(self._algo_refine.refinement_strategy(),self._gmesh)
        #print ( "geom_bfun_set._nbfuns_active = %d\n" % geom_bfun_set._nbfuns_active )
        #print ( "geom_bfun_set._nbfuns_total = %d\n" % geom_bfun_set._nbfuns_total )
        
        #nbfuns = len(geom_bfun_set._bfuns)
        #print("nbfuns=%d\n"%nbfuns)
        #for i in range(0,nbfuns):
        #    print(geom_bfun_set._bfuns[i])
        #print("nbfuns=%d\n"%nbfuns)
        #print ( geom_bfun_set._bfuns[0]._gcells )
        #for i in range(0,geom_bfun_set._nbfuns_active): #self._nbfuns_total
        #    print("fen_label = %d, no. of gcells = %d"%(geom_bfun_set._bfuns[i]._fen._id,len(geom_bfun_set._bfuns[i]._gcells)))
        
        #    #print("%d, nbfuns = %d\n"%(i,len(bfuns.)))
        #print(geom_bfun_set._bfun_dofparam_ids)
        
        self._geometry = FIELD_VECTOR("geometry",geom_bfun_set,zero)
        
        raise NotImplementedError
    
    def solve(self):
        pass
    
    def geometry(self):
        return self._geometry
    
    def u(self):
        return self._u
    
    def gmesh(self):
        return self._gmesh
    
    def adapt(self):
        pass
    