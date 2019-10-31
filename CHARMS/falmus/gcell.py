from abc import ABC, abstractmethod
from famuls import CHECK
from conn import CONN_SOLID_8

class GCELL(ABC):
    
    def __init__(self):        
          self._gcell_group = None
          self._id = None

    @abstractmethod
    def type_name(self):
        pass    
    
    @staticmethod
    def make_gcell(gcell_type, implementation, fens ):
        if ( "solid_h8" in gcell_type ):
            return GCELL_SOLID_H8(fens)
        else:
            raise Exception('make_gcell Fails - Element Type Not Defined\n')
    
    # Return the gcell group to which this gcell belongs.
    # A gcell may belong to a single group only.     
    def gcell_group(self):
        return self._gcell_group
            
    def set_gcell_group(self,gcell_group):
        self._gcell_group = gcell_group
      
    @abstractmethod
    def manifold_dim(self):
        pass
    
    # Return a pointer to the connectivity this gcell is based upon.
    # Each gcell is based upon some connectivity (but connectivity
    # does not necessarily have to be associated with a gcell: consider
    # the edges of a triangular gcell, they are not gcells, but they
    # are connectivities).  Each gcell specialization uses a particular
    # connectivity to express its manifold dimension, number of connected
    # nodes, and the refinement pattern.
    @abstractmethod
    def conn(self):
        pass

    # Evaluate the basis function set associated with the evaluation point.
    #  It is assumed that the gcell on which this methods gets invoked
    #  is at the top of the hierarchy.  The reason is that the basis
    #  function evaluation proceeds from the finest level (top) to the
    #  coarsest level (bottom).
    #  The standard finite element basis functions defined
    #  over the gcell specialization (for instance, the tri-linear functions
    #  over an 8-node hexahedron) need to be evaluated for each cell in the
    #  hierarchy.  The basis functions (and their derivatives) are
    #  computed in the natural coordinate system.  When the whole hierarchy
    #  had been traversed, the Jacobian transformation is evaluated
    #  and the derivatives are transformed into the given (Cartesian)
    #  coordinate system. 
    @abstractmethod
    def eval_bfun_set(self,evalpt):
        pass
    
    
    # ==========================================================
    #  methods implement the gcell refinement hierarchy 
    # ==========================================================
    
    # At which level in the hierarchy is the cell?
    # The coarsest level is zero (0).
    def level(self):
        l = 0
        gcell = self        
        while ( gcell.parent() is not None ):
            
            l += 1 
            gcell = gcell.parent()
  
        return l
  
    # Return the parent of this gcell.  May be null if the
    # cell has no parent (meaning it is at the bottom).
    @abstractmethod
    def parent(self):
        pass
    
    
    # Get the i-th child of the cell.
    @abstractmethod
    def child(self,i):
        pass
    
    # Get the number of children.  If the gcell is not refined, zero should be
    # returned; otherwise the number of child gcells is returned
    # (depending on the type of the gcell).
    @abstractmethod
    def nchildren(self):
        pass

    # Has the gcell been divided into children?  Note we're not asking if they
    # *could* exist, rather whether they actually exist.
    @abstractmethod
    def divided(self):
        pass
  
    # Physically divide the gcell by creating its children (and all the nodes
    # these children connect on the higher level).  Needs to be defined
    # by the implementation of the specialization.
    @abstractmethod
    def divide(self,ref_ctx):
        pass
    
    # Return the set of detail functions that are needed for the refinement
    #  of the function associated with the node on input.
    #  This function may be called *ONLY* if the gcell had been divided
    #  into children before.
    @abstractmethod    
    def detail_set(self,fen_of_bfun_to_refine,rf):
        pass
  
    # Return the set of all refinement functions (i.e. detail functions plus
    #  the one private function) that are needed for the refinement
    #  of the function associated with the node on input.
    #  This function may be called *ONLY* if the gcell had been divided
    #  into children before.  
    @abstractmethod    
    def complete_refinement_set(self,fen_of_bfun_to_refine,rf):
        pass

    # Map a parametric point xi, eta, theta (one, two, or all three of these
    # coordinates may be used) into a child of this gcell.
    @abstractmethod    
    def map_to_child(self,param_loc,child,child_param_loc):
        pass
    
  
    # Given a parametric location, return the parent and the parametric coordinates
    #  to which this location maps in the parent of the this cell.
    #  If the cell has no parent, false is returned and the output parameters are undefined.
    #  If successful, *parent is the parent of the cell, and
    #  true is returned.
    @abstractmethod    
    def map_to_parent(self,param_loc, parent, parent_param_loc):
        pass


    # Map a finite element node to the parametric coordinate in the
    #  gcell given as argument.  If the finite element node happens
    #  not to be referenced by the connectivity of the cells, an
    #  exception EXCEPTION_BAD_VALUE is thrown.
    @abstractmethod    
    def map_fen (fen, param_loc):
        pass


    # Is this the topmost gcell that is used by any basis function in the field?
    def is_topmost_gcell_in_field(self,field):
        pass
    
    
    # Find the topmost gcell that is used by any basis function in the
    #  field given as argument.  If no such gcell exists, null (0) is returned.
    #  If non-null is returned, *txi, *teta, and *ttheta return the parametric
    #  coordinates of the point xi, eta, theta mapped into the topmost gcell.
    #  Note: This is a generic procedure that relies on the implementations
    #  of map_to_child() and map_to_parent() in the specializations of the gcell.
    def map_to_topmost_gcell_in_field(self, field, param_loc, mapped_to_param_loc):
        pass

  
    # Find the root gcell of the refinement hierarchy whose member `this' is.
    # (It's the Adam or Eve ;)  
    def root(self):
        
        r = self        
        while (r.parent()):
            r = r.parent()
    
        return r

    # =======================
    # debugging methods 
    # =======================
  
    def set_id(self,id):
        CHECK (self._id == None) #,EXCEPTION_ILLEGAL_USE,;);
        self._id = id
  
    def id(self):
        return self._id
    

    
class GCELL_GROUP:
    
    class gcell_enumerator_t():        
        
        def __init__(self,gcells):
            self._gcells = gcells
            
        def __iter__(self):
            for gcell in self._gcells:
                yield gcell
            
    def gcell_enumerator(self):
        e = self.gcell_enumerator_t(self._gcells)
        return e
    
            
    def __init__(self,name,db,gsubmesh):
        self._name       = name
        self._gsubmesh   = gsubmesh
        self._gcell_type = ""
        self._gcells = []
        
        # type
        path = "algorithms/gcell_groups/" + self._name + "/type"
        self._gcell_type = db.DB_GET_STRING(path)
        #print(self._gcell_type)
        
        # conn_file
        path = "algorithms/gcell_groups/" + self._name + "/conn_file"
        gcell_file = db.DB_GET_STRING (path)
        CHECK(self.read_gcell_file(gcell_file))#, EXCEPTION_BAD_VALUE ,;);


    def add(self,gcell):        
        CHECK (self._gcell_type == gcell.type_name()) #, EXCEPTION_BAD_VALUE,;);
        self._gcells.append(gcell)
        gcell.set_gcell_group(self)


    def read_gcell_file (self,file):
        
        read_gcell_flag = False
        
        content=[]
        with open(file,'r') as f:
            for line in f:
                content.append(line)
            
        for eid,line in enumerate(content):
            line = line.strip().split()
            
            fens = []
            for id in line:
                id = int(id)
                fen = self._gsubmesh.gmesh().find_fen(id)
                if fen is None:
                    print("Node " + str(id) + " not found\n")
                    CHECK (fen is not None) #, EXCEPTION_NULL_PTR,;);
                fens.append(fen)
        
            gcell = GCELL.make_gcell(self._gcell_type, "default", fens)
            if (gcell is None):
                print("While reading \"" + file + "\": making " + str(self._gcell_type))
                CHECK (gcell is not None)#, EXCEPTION_NULL_PTR,;);
            
            gcell.set_id(eid+1)
            self.add(gcell)
    
        read_gcell_flag = True
        return read_gcell_flag 

# ---------------------------------------------------------------------------


class GCELL_SOLID_H8(GCELL):
    
    NCHILDREN   = 8
    NDETAILFENS = 19 # mid-edge nodes, mid-face nodes, mid-volume node
    NEDGES      = 12
    NFACES      = 6
    _child_map  = [ [ [-1,-1,-1], [0,0,0] ],
                  [ [ 0,-1,-1], [1,0,0] ],
                  [ [ 0, 0,-1], [1,1,0] ],
                  [ [-1, 0,-1], [0,1,0] ],
                  [ [-1,-1, 0], [0,0,1] ],
                  [ [ 0,-1, 0], [1,0,1] ],
                  [ [ 0, 0, 0], [1,1,1] ],
                  [ [-1, 0, 0], [0,1,1] ] ]
    
    # param coords correponding to vertices
    g3dhex_vertex_param_coord = [ [ -1, -1, -1 ],
                                  [  1, -1, -1 ],
                                  [  1,  1, -1 ],
                                  [ -1,  1, -1 ],
                                  [ -1, -1,  1 ],
                                  [  1, -1,  1 ],
                                  [  1,  1,  1 ],
                                  [ -1,  1,  1 ] ]


    
    def __init__(self,fens):   
        GCELL.__init__(self)
        self._conn   = CONN_SOLID_8()
        
        self._conn._fens = fens
        
        CHECK( len(fens) == self._conn.nfens() )
        self.TYPE_NAME = "solid_h8"
        
        self._parent = None
        self._child  = []
    
    def manifold_dim(self):
        return self._conn.manifold_dim()
        
    def type_name(self):
        return self.TYPE_NAME


    # Return connectivity.
    def conn(self):
        return self._conn

    def eval_bfun_set (self,evalpt):
        raise NotImplementedError

    def map_fen (self,fen, param_loc):
        for j in range(0,self._conn.nfens()):
            if ( fen.id() == self._conn.fen(j).id() ):
                param_loc[0] = self.g3dhex_vertex_param_coord[j][0]
                param_loc[1] = self.g3dhex_vertex_param_coord[j][1]
                param_loc[2] = self.g3dhex_vertex_param_coord[j][2]
                
                return
        
        raise Exception  # the node is unknown to this cell; raise hell

    def map_to_child(self, param_loc, refchild, child_param_loc):
        
        CHECK( isinstance(refchild,list) == True )
        # note: 
        #  only child_param_loc is modified
        #  refchild is NOT modified !!!!!!!!!!!!!!!
        # but refchild[0] is
        
        
        if (self.nchildren() == 0):
            return False

        xi    = param_loc(0)
        eta   = param_loc(1)
        theta = param_loc(2)        
        
        if (xi < 0):
            
            if (eta < 0):
                if (theta < 0):
                    refchild[0] = self._child[0]
                    child_param_loc[0] = 2 * xi    + 1
                    child_param_loc[1] = 2 * eta   + 1
                    child_param_loc[2] = 2 * theta + 1
                else:
                    refchild[0] = self._child[4]
                    child_param_loc[0] = 2 * xi    + 1
                    child_param_loc[1] = 2 * eta   + 1
                    child_param_loc[2] = 2 * theta - 1
            
            else: # if eta > 0
                
                if (theta < 0):
                    refchild[0] = self._child[3]
                    child_param_loc[0] = 2 * xi    + 1
                    child_param_loc[1] = 2 * eta   - 1
                    child_param_loc[2] = 2 * theta + 1
                else:
                    refchild[0] = self._child[7]
                    child_param_loc[0] = 2 * xi    + 1
                    child_param_loc[1] = 2 * eta   - 1
                    child_param_loc[2] = 2 * theta - 1
        
        else: # if xi > 0
            
            if (eta < 0):
                
                if (theta < 0):
                    
                    refchild[0] = self._child[1]
                    child_param_loc[0] = 2 * xi    - 1
                    child_param_loc[1] = 2 * eta   + 1
                    child_param_loc[2] = 2 * theta + 1
                
                else:
                    
                    refchild[0] = self._child[5]
                    child_param_loc[0] = 2 * xi    - 1
                    child_param_loc[1] = 2 * eta   + 1
                    child_param_loc[2] = 2 * theta - 1
      
            else: # if eta > 0
                
                if (theta < 0):
                    refchild[0] = self._child[2]
                    child_param_loc[0] = 2 * xi    - 1
                    child_param_loc[1] = 2 * eta   - 1
                    child_param_loc[2] = 2 * theta + 1
                else:
                    
                    refchild[0] = self._child[6]
                    child_param_loc[0] = 2 * xi    - 1
                    child_param_loc[1] = 2 * eta   - 1
                    child_param_loc[2] = 2 * theta - 1
  
        
        return True
    
    
    def map_to_parent(self,param_loc, refparent, parent_param_loc):
        
        CHECK( isinstance(refparent,list) == True )
        
        if ( self._parent is not None ):
            
            xi = param_loc(0)
            eta = param_loc(1)
            theta = param_loc(2)
            
            for j in range(0,GCELL_SOLID_H8.NCHILDREN):
                
                if ( self == self._parent.child(j)):                
                    LO = 0
                    HI = 1
                    parent_param_loc[0] = (1 -    xi) / 2 * self._child_map[j][LO][0] + (1 +    xi) / 2 * self._child_map[j][HI][0]
                    parent_param_loc[1] = (1 -   eta) / 2 * self._child_map[j][LO][1] + (1 +   eta) / 2 * self._child_map[j][HI][1]
                    parent_param_loc[2] = (1 - theta) / 2 * self._child_map[j][LO][2] + (1 + theta) / 2 * self._child_map[j][HI][2]
                    refparent[0] = self._parent
                    return True
            
            # end loop
            
            return False # found no parent

    def map_to_parent_fen (self,fen):
        raise NotImplementedError # not implemented in famuls
        
    def divide (self,ref_ctx):
        raise NotImplementedError
        
    def detail_set(self,fen_of_bfun_to_refine, rf):
        raise NotImplementedError
                
    def complete_refinement_set(self,fen_of_bfun_to_refine, rf):
        raise NotImplementedError
        
    def parent(self):
        return self._parent
    
    def child(self,i):
        return self._child[i]
    
    def nchildren(self):
        if ( ( len(self._child) == self.NCHILDREN ) and ( self._child[0] is not None ) ): # We're not checking every child: if there is one, there have to be all of them
            return self.NCHILDREN 
        else:
            return 0
    
    def divided(self):
        return (self.nchildren() > 0)    
    
    def eval_ns(self, child, evalpt, xi, eta, theta, dxiparent_dxichild):
        raise NotImplementedError
                  
