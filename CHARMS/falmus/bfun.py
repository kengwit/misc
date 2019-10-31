from famuls import CHECK
from abc import ABC, abstractmethod
from gmesh import GSUBMESH,GMESH 
from fen import FEN
from enum import Enum
from bfun_dofparam_pair_id import INVALID_BFUN_DOFPARAM_PAIR_ID

class BFUN_SET:
    
    class REFINEMENT_STRATEGY(Enum):    
        QUASI_HIERARCHICAL = 0
        TRUE_HIERARCHICAL  = 1
    
    USE_PU = False
    
    class gsubmesh_enumerator_t():        
        
        def __init__(self,gsubmeshes):
            self._gsubmeshes = gsubmeshes
            
        def __iter__(self):
            for gsubmesh in self._gsubmeshes:
                yield gsubmesh
            
    
    # Needed to ensure that iteration proceeds over active functions only.
    class bfun_enumerator_t():        
        
        def __init__(self,bfuns):
            self._bfuns = bfuns
            
        def __iter__(self):
            for bfun in self._bfuns:
                if ( ( bfun is not None ) and ( bfun.is_active() ) ):
                    yield bfun 
    
    class bfun_enumerator_all_t():        
        
        def __init__(self,bfuns):
            self._bfuns = bfuns
            
        def __iter__(self):
            for bfun in self._bfuns:
                yield bfun 

    # Enumerator of submeshes.
    def gsubmesh_enumerator(self):
        e = self.gsubmesh_enumerator_t(self._gsubmeshes)
        return e

    # Enumerator of active basis functions.
    def bfun_enumerator(self):
        #print("self._bfuns")
        #print(self._bfuns[0] is not None) # this gives True
        #print(not self._bfuns[0].is_active() ) # this gives False !!!!!!!!
              
        e = self.bfun_enumerator_t(self._bfuns)
        return e
    
    def bfun_enumerator_all(self):
        e = self.bfun_enumerator_all_t(self._bfuns)
        return e
    
    def __init__(self,refinement_strategy,*args):
        
        if ( len(args) == 0 ):
            #print('zero args')
            self._in_constructor = True
            self._refinement_strategy = refinement_strategy
            self._gmesh      = None
            self._gsubmeshes = []
            self._bfuns      = []
            self._nbfuns_active     = 0
            self._nbfuns_total      = 0
            self._bfun_dofparam_ids = []
            self._in_constructor    = False
            self._mem_raise_size = 10
            
            if BFUN_SET.USE_PU:
                self._is_pu = False
        
        elif isinstance( args[0], GSUBMESH ):
            #print('gsubmesh args')
            
            gsubmeshes = args[0]
            
            self._mem_raise_size = 10
            self._in_constructor = True
            self._refinement_strategy = refinement_strategy
            
            self._gsubmeshes = []
            for gsm in gsubmeshes:
                CHECK (gsm is not None) #, EXCEPTION_NULL_PTR,;);
                self._gsubmeshes.append(gsm)
                
            self._gmesh = self._gsubmeshes[0].gmesh() # `local' field: constructed from submeshes
            self._bfuns         = []
            self._nbfuns_active = 0
            self._nbfuns_total  = 0
            self._bfun_dofparam_ids = []
            
            igcells = self.build_gcell_list()
            CHECK (len(igcells) != 0 ) #, EXCEPTION_BAD_VALUE,;);
            self.build_from_gcells (igcells, 0)
            
            if BFUN_SET.USE_PU:
                self._is_pu = False
            
            self._in_constructor = False
            
            
        elif isinstance( args[0], GMESH ):
            #print('gmesh args')
            
            gmesh = args[0]
            
            self._mem_raise_size = 10
            self._in_constructor = True
            self._refinement_strategy = refinement_strategy
                        
            self._gmesh = gmesh #  // `global' field: constructed from the mesh
            self._gsubmeshes = []
            
            for gsm in gmesh._gsubmeshes:
                self._gsubmeshes.append(gsm)
                
            self._bfuns         = []
            self._nbfuns_active = 0
            self._nbfuns_total  = 0
            self._bfun_dofparam_ids = []
            
            igcells = self.build_gcell_list()
            CHECK (len(igcells) != 0 ) #, EXCEPTION_BAD_VALUE,;);
            self.build_from_gcells(igcells, 0)
            
            if BFUN_SET.USE_PU:
                self._is_pu = False
            
            self._in_constructor = False
            
        else:
            raise NotImplementedError
            
        # end if

    # end def 
    
    @staticmethod
    def make(refinement_strategy,*args):
        sh = BFUN_SET(refinement_strategy,*args)

        if BFUN_SET.USE_PU:
            sh.make_pu()

        return sh
            
    # Return one if refinement strategy is true hierarchical. Otherwise it returns 0.
    def refinement_strategy(self):
        return self._refinement_strategy
        
    
    # Returns true if basis function passed as argument exists as inactive basis function in 
    # basis function set.
    def in_inactive_set (self,bfun):
        return ( (self.dofparam_id(bfun.fen()) >= 0) and ( not self.in_active_set(bfun) ) )
  
  
    # Returns true if basis function passed as argument exists as active basis function in 
    # basis function set.
    def in_active_set (self,bfun):
        dpid = self.dofparam_id(bfun.fen());
        return ( ( dpid >= 0) and (dpid < self._nbfuns_active ) )
  
    # Is this a partition-of-unity basis function set?
    def is_pu(self):
        if BFUN_SET.USE_PU:
            return self._is_pu 
        else:
            return False
        
  
    def make_pu(self):
        if BFUN_SET.USE_PU:
            raise NotImplementedError  
        else:
            pass
        
  
    # Get the bfun/dofparam id.  It is the (opaque) id that fields
    # use to access pairs.  All fields based on the same bfun_set
    # will have the same order of pairs, and hence all may use the
    # dofparam identifier generated by the bfun set. INVALID_DOFPARAM_ID
    # is returned if basis function is not active.
    #
    # Remark : Bfun has to be associated with given basis function set.
    def dofparam_id(self,*args):
        
        #print("in dofparam_id, args is: ", end='')
        #print(args,end='')
        if isinstance( args[0], BFUN ):
            #print(". It is an BFUN, fen label = ",end='')
            bfun = args[0]
            #print(bfun._fen.id())
            CHECK (bfun.bfun_set() == self) #, EXCEPTION_BAD_ACCESS,;);
            
            if ( bfun.is_active() is not True ):
                return INVALID_BFUN_DOFPARAM_PAIR_ID # RAW
            else:
                return self.dofparam_id( bfun.fen() )
        
        elif isinstance( args[0], FEN ):
            #print(". It is an FEN, fen label = ",end='')
            fen = args[0]
            #print(fen.id())
            if ( fen.uniqobjid() >= len( self._bfun_dofparam_ids ) ):
                #print("fen.uniqobjid() >= len( self._bfun_dofparam_ids )\n" )
                #raise Exception
                #return sys.maxsize 
                return INVALID_BFUN_DOFPARAM_PAIR_ID
            else:
                return self._bfun_dofparam_ids[ fen.uniqobjid() ]
            
        else:
            print("Error, args is: ",end='')
            print(args)
            raise NotImplementedError            
            
    
    def is_present(self,bfun):
        return (self.dofparam_id( bfun.fen() ) != INVALID_BFUN_DOFPARAM_PAIR_ID)
        
    
    def insert_bfun (self,bfun):
        
        # be careful here, it is not simply appending
        if ( len( self._bfuns ) <= self._nbfuns_total ):
            start = len( self._bfuns )
            end   = start + self._mem_raise_size
            for k in range(start, end):
                self._bfuns.append(None)
  
        self._bfuns[ self._nbfuns_total ] = bfun
        self._nbfuns_total += 1
                
        
        if ( len( self._bfun_dofparam_ids ) <= bfun.fen().uniqobjid() ):
            
            p_size = len(self._bfun_dofparam_ids) # size before resizing
            
            new_total_size = bfun.fen().uniqobjid()+self._mem_raise_size
            while True: # infinite loop; resize until you gete the size you want
                self._bfun_dofparam_ids.append(None)
                if (len(self._bfun_dofparam_ids)==new_total_size):
                    break
                        
            for k in range(p_size,len(self._bfun_dofparam_ids)):
                self._bfun_dofparam_ids[k]= INVALID_BFUN_DOFPARAM_PAIR_ID
        
        self._bfun_dofparam_ids[bfun.fen().uniqobjid()] = self._nbfuns_total-1
        
        
        
        
        
    def swap(self,temp_loc1, temp_loc2):
        temp_bfun   = self._bfuns[temp_loc1]
        self._bfuns[temp_loc1] = self._bfuns[temp_loc2]
        self._bfuns[temp_loc2] = temp_bfun
        
        temp_dpid = self._bfun_dofparam_ids[ self._bfuns[temp_loc1].fen().uniqobjid() ]
        
        self._bfun_dofparam_ids[ self._bfuns[temp_loc1 ].fen().uniqobjid()] = \
            self._bfun_dofparam_ids[self._bfuns[temp_loc2].fen().uniqobjid()] 
            
        self._bfun_dofparam_ids[self._bfuns[temp_loc2].fen().uniqobjid()] = temp_dpid
          
  
    def build_from_gcells (self, igcells, bfun_set):
        # here, we basically get every cell that is connected to a node
        
        #inactive_bfun = [] # not used
        fen_to_gcells_map = {} # map of a list, i.e. fen --> [gcell1, ...]
                
        for gcell in igcells:
            conn = gcell.conn()
            nfens = conn.nfens()
            
            for j in range(0,nfens):
                fen = conn.fen(j)
                
                if ( fen not in fen_to_gcells_map ):
                    fen_to_gcells_map[ fen ] = [] # empty list
                # end if
                
                fen_to_gcells_map[ fen ].append( gcell )
        
        # now iterate the map and set up the basis functions
        for fen in fen_to_gcells_map:
            # Decide whether that node is to be included:
            gcells = fen_to_gcells_map[ fen ]
            bfun = BFUN_FE (fen, gcells, self ) # RAW : should be replaced with a make function.
            #print("setting bfun active = True")
            bfun.set_is_active(True)
            self.add(bfun)
      
        
        
    
    def build_gcell_list(self):
        igcells = []
        #print(self._gsubmeshes)
        e0 = self.gsubmesh_enumerator()
        #for gsubmesh in self._gsubmeshes:
        for gsubmesh in e0:
            
            e = gsubmesh.gcell_group_enumerator()  
            #for gcell_group in gsubmesh._gcell_groups:
            for gcell_group in e:
                
                e1 = gcell_group.gcell_enumerator()      
                #for gcell in gcell_group._gcells: # all elements in the group
                for gcell in e1: # all elements in the group
                    igcells.append(gcell)
    
        return igcells
    
    def build_bfun_dofparam_ids(self):
        raise NotImplementedError

    def nbfuns(self):
        return self._nbfuns_active
  
    # add a basis function to basis function set. This function can be used to change the attribute 
    # of a basis function from active to inactive or vice versa.  
    def add (self,bfun):
        CHECK (self._in_constructor == True) #,EXCEPTION_ILLEGAL_USE,;); //RAW : now we are adding bfun in clone 
        CHECK(bfun.bfun_set() == self)#, EXCEPTION_ILLEGAL_USE,;);
        #print('In ADD\n')
        if (not self.is_present(bfun)):
            
            self.insert_bfun(bfun)
            
            if ( bfun.is_active() ): # if basis function has to be  active. 
     
               if ( not self.in_active_set(bfun) ): # but it was grouped among inactive
                    self.swap( self.dofparam_id( bfun.fen() ),self._nbfuns_active)
                    self._nbfuns_active += 1


            else: # if a basis function has to be inactive 
                
                if ( self.in_active_set(bfun) ): # but it was grouped among actives 
                   self.swap( self.dofparam_id(bfun.fen()), self._nbfuns_active-1)
                   self._nbfuns_active -= 1
            
            # endif
            
        # end if
        

  
# -------------------------------------------------------------------------
class BFUN(ABC):
    

    
    def __init__(self,fen,bfun_set,gcells):
        self._gcells = []
        for gcell in gcells:
            self._gcells.append(gcell)
            
        self._is_active = False
        self._bfun_set  = bfun_set
        self._fen       = fen
        self._ref_set   = []
        self._beta = 1 #non-PU
  
    # Get the finite element node associated with this basis function.
    def fen(self):
        return self._fen
    
    # Get the bfun_set associated with current basis function.
    # Note : Correspondance between basis function set and basis function is one to one. 
    def bfun_set(self):
        return self._bfun_set
    
    # Clone a basis function.
    @abstractmethod
    def clone(self, bfun_set):
        pass

    # Return true if basis function is active.
    def is_active(self):
        return self._is_active
        
    # set the flag that says wether a function is active or not. 
    def set_is_active(self,is_active):
        self._is_active = is_active        
        
    # Return number of gcells covered (at least partially) by the bfun.  
    def ngcells(self,):
        return len( self._gcells )
             
    # Return a gcell.
    def gcell(self,j):
        return self._gcells[j]



    #@ TODO this needs to be refined
    # return the detail_set or refinement set  of basis function depending on the refinement strtegy of 
    # basis function set.  
    def refinement_set(self,s):
        raise NotImplementedError
        if ( self.bfun_set().refinement_strategy() == self.bfun_set().REFINEMENT_STRATEGY.TRUE_HIERARCHICAL ):
            self.detail_set(s)
        else:
            self.complete_refinement_set(s)
            
        return True

    #@ TODO this needs to be refined    
    def complete_refinement_set(self,s):
        raise NotImplementedError
        s = []
        for gc in self._gcells:
            if ( gc.divided() ):
                gc.complete_refinement_set(self.fen(),s)  
        return

    #@ TODO this needs to be refined    
    def detail_set(self,s):
        raise NotImplementedError
        s = []
        for gc in self._gcells:
            if ( gc.divided() ):
                gc.complete_refinement_set(self.fen(),s)  
        return
      
    # Return true if basis function was refined.
    def is_refined(self):
        #for j in range(0,self.ngcells()):
        #    if ( not self.gcell(j).divided() ):
        #        return False
            
        raise NotImplementedError
            
        
# -------------------------------------------------------------------------
class BFUN_FE(BFUN):
    
    def __init__(self,fen, gcells, bfun_set):
        BFUN.__init__(self,fen,bfun_set,gcells)
        
        
        
    def clone (bfun_set: BFUN_SET):
        raise NotImplementedError

        
#    def __init__(self,refinement_strategy,gsubmeshes,gmesh):
#        
#        if ( gmesh is None and gsubmeshes is not None):            
#            print('BFUN_SET1')
#        elif ( gmesh is not None and gsubmeshes is None):            
#            print('BFUN_SET2')
#        
#        else:
#            print('BFUN_SET3')
#
#            self._in_constructor = True
#            self._refinement_strategy = refinement_strategy
#            self._gmesh = None
#            self._gsubmeshes = []
#            self._bfuns = []
#            self._nbfuns_active = 0
#            self._nbfuns_total = 0
#            self._bfun_dofparam_ids = []
#            self._in_constructor = False
#            self._mem_raise_size = 10
#            #if defined(USE_PU) && USE_PU
#            self._is_pu = False
#            #endif  
        
    

        
    #@staticmethod
    #def make(refinement_strategy,gsubmeshes,gmesh):
#        if ( gmesh is None and gsubmeshes is not None):            
#            sh = BFUN_SET(refinement_strategy,gsubmeshes,None)
#            #sh.make_pu()
#            return sh
#        elif ( gmesh is not None and gsubmeshes is None):            
#            sh = BFUN_SET(refinement_strategy,None,gmesh)
#            #sh.make_pu()
#            return sh
#        
#        else:
#            sh = BFUN_SET(refinement_strategy,None,None)