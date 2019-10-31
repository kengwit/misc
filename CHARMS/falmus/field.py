from bfun import BFUN_SET
from abc import ABC, abstractmethod
from famuls import CHECK
#from typing import TypeVar, Generic, overload
from bfun_dofparam_pair_id import INVALID_BFUN_DOFPARAM_PAIR_ID


class FIELD_BASE(ABC):
    
    
        
#    def __init__(self,name,bfun_set):
#            self._name     = name
#            self._bfun_set = bfun_set

    def __init__(self,*args):            
        #print(args)
        if ( len(args) == 0 ):
            self._name = ""
            self._bfun_set = None
            self.NCOMPONENTS = 0
        
        elif ( isinstance( args[0], str ) and isinstance( args[1], BFUN_SET ) and isinstance( args[2], int ) ):
            
            self._name     = args[0]
            self._bfun_set = args[1]            
            self.NCOMPONENTS = args[2]
            
            
        else:
            
            raise NotImplementedError

    # Provide the caller with an opaque identifier for a degree of freedom
    #  associated in the field with given node.  If this identifier is passed
    #  to *any* field based on the same basis functions set, the field will be
    #  able to return the value of the degree of freedom in constant time (i.e. quickly). 
    @abstractmethod
    def dofparam_id(self,fen):
        pass
    
    # Return the basis function set on which this field is based.
    def bfun_set(self):
        return self._bfun_set
  
    # Get the pairs active (ie. non-zero basis function) over a gcell.
    @abstractmethod
    def pairs_active_over_gcell(self,gcell):
        pass
  
    # Is the field active over cell?
    @abstractmethod
    def active_over(self,gcell):
        pass
    
    
    # Return the bfun associated with the identifier on input.   This operation will be done
    # efficiently (in constant time).
    @abstractmethod
    def bfun(self,dofparam_id):
        pass
  
    # Return the name.
    def name(self):
        return self._name
            
class FIELD_PAIR:
    
    NUM_COMPONENTS = 0
    
    def __init__(self, bfun, initial_value):
        self._bfun        = bfun
        self._dofparam    = initial_value
        self._constrained = 0
        
    def bfun(self):
        return self._bfun
    

class FIELD( FIELD_BASE ):

    def __init__(self,name,bfun_set,initial_value,*args):
        
        if ( len (args) == 1 and isinstance( args[0], int ) ):
            ncomps = args[0]
            #print("ncomps=%d\n"%ncomps)
            FIELD_BASE.__init__(self,name,bfun_set,ncomps)    
        else:
            print("must have NCOMPONENTS input")
            raise NotImplementedError
        
                
        self._gcell_to_bfun_dofparam_id_map = {}
        self._pairs = [None] * bfun_set.nbfuns()        
        
        # build the list of pairs - use only active basis functions !!!!!!!!
        e = bfun_set.bfun_enumerator()      
        #print("e=",end='')
        #print(e)
        
        # " BELOW LOOP IS NOT EXCECUTED !!!!!!!!!!!! "
        
        for bfun in e: 
            #print("--------------------------------------------------\n")
            #print("counter = %d, in loop create field pair"%counter)
            field_pair = FIELD_PAIR(bfun,initial_value)
            field_pair.NUM_COMPONENTS = self.NCOMPONENTS
            #print("field_pair.NUM_COMPONENTS = %d\n"%field_pair.NUM_COMPONENTS)
            #print("done create field pair")
            
            dpid = bfun_set.dofparam_id(bfun)
            
            #print("dpid = ", end='')
            #print(dpid,end='')
            #print(", field_pair = ",end='')
            self._pairs[ dpid ] = field_pair
                    
        
        #print("done")
        # and set up the map from gcells to active (ie. non-zero bfun) field pairs 
        self.build_gcell_to_bfun_dofparam_id_map()
        
    
    def build_gcell_to_bfun_dofparam_id_map(self):
        _M = self._gcell_to_bfun_dofparam_id_map
        _M = {}
        # Add empty lists to the map for all gcells active in the field
        #print(self._pairs)
        
        for j in range(0,len(self._pairs)):
            bfun = self._pairs[j].bfun()
            for k in range(0,bfun.ngcells()):
                if ( bfun.gcell(k) not in _M ):
                     _M[ bfun.gcell(k) ] = set() # empty set
            
        # Now for each pair ...
        for j in range(0,len(self._pairs)):
            pair = self._pairs[j]
            bfun = pair.bfun()
            for k in range(0,bfun.ngcells()): # ... loop over its cells            
                gcell = bfun.gcell(k)                                
                CHECK (gcell in _M) #, EXCEPTION_BAD_ACCESS,;);
                
                # ------------------------------------
                # First propagate your influence to 
                # your ancestors
                # ------------------------------------
                parent = gcell.parent()
                
                while ( parent is not None ):
                    
                    if ( parent in _M ):
                        _M[ parent ].add( j )
                        
                    parent = parent.parent()
                
                # end while        
                
                # ------------------------------------
                # Now for yourself and your children
                # ------------------------------------
                s = []
                s.append(gcell) # push
                while ( len(s) != 0 ):
                    
                    agcell = s[-1] # top element                        
                    
                    s.pop() #  working on this gcell
                    
                    if ( agcell in _M ): # This cell is active in the field
                        _M[ agcell ].add( j )
                        
                    # end if
                    
                    # work on the children too
                    for l in range(0,agcell.nchildren()): 
                        
                        s.append(agcell.child(l))
                
                # end while
                
            # end loop over each gcell in pair support
        # end loop over each pair
        
        
    def pairs_active_over_gcell(self,gcell):
        raise NotImplementedError
        
        if ( gcell in self._gcell_to_bfun_dofparam_id_map ):
            bfun_dofparam_id_set = self._gcell_to_bfun_dofparam_id_map[gcell]
            return bfun_dofparam_id_set
        else:
            return set()
            
        
    def active_over(self,gcell):
        raise NotImplementedError
        
        if ( gcell in self._gcell_to_bfun_dofparam_id_map ):
            bfun_dofparam_id_set = self._gcell_to_bfun_dofparam_id_map[gcell]
            return ( len(bfun_dofparam_id_set) > 0 )
        else:
            raise False
            
      
    def bfun(self,dofparam_id):  
        CHECK (dofparam_id != INVALID_BFUN_DOFPARAM_PAIR_ID ) #, EXCEPTION_BAD_ACCESS,;);
        return self._pairs[dofparam_id].bfun()
  

    def dofparam_id(self,fen):
        bfun = self.bfun_set().fen_to_bfun(fen)
        if (bfun is None):
            return INVALID_BFUN_DOFPARAM_PAIR_ID
        else:
            return self.bfun_set().dofparam_id(bfun)
  

def FIELD_SCALAR(name,bfun_set,initial_value):
    
    fs = FIELD(name,bfun_set,initial_value,1)
    
    return fs

def FIELD_VECTOR(name,bfun_set,initial_value):
    
    fv = FIELD(name,bfun_set,initial_value,3)
    
    return fv

    
    
    
