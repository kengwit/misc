from famuls import CHECK
from abc import ABC, abstractmethod
from enum import Enum
from typing import TypeVar, Generic

class CONN_MANIFOLD_DIM(Enum):
  CONN_0_MANIFOLD = 0
  CONN_1_MANIFOLD = 1
  CONN_2_MANIFOLD = 2
  CONN_3_MANIFOLD = 3
        


class CONN_BASE(ABC):
    
    def __init__(self):
        pass
    
    # Return the number of nodes that this connectivity binds.
    @abstractmethod
    def nfens(self):
        pass
    
    # Return the i-th node.
    # Return the i-th node for modification.
    @abstractmethod
    def fen(self):
        pass
    
    # What is the manifold dimension of the conn?
    # (Point=0, line=1, surface=2, or solid=3.) 
    @abstractmethod
    def manifold_dim(self):
        pass

  
    # Get connectivity code.
    @abstractmethod
    def conn_code (self):
        pass

  
    # Return the number of refinement nodes that this connectivity uses
    # in a refinement step.
    @abstractmethod
    def nreffens (self):
        pass
  
    # Get a refinement node.
    # The refinement nodes are oriented for the connectivity on
    # input.  Therefore, first try to see if there is any topological
    # transformation that may be applied to the connectivity on input
    # that would lead to a successful match with self.  If that is
    # the case, the same topological operations need to be applied
    # to the index into vector of refinement nodes indx, and
    # ref_fens[transf(indx)] is returned; otherwise, that is if match
    # is not possible, null is returned.
    @abstractmethod
    def get_ref_fen (conn, ref_fens, indx):
        pass
        
    @abstractmethod
    def clone ():
        pass
        
        
MANIFOLD_DIM = TypeVar('CONN_MANIFOLD_DIM')
NFENS        = TypeVar('Int')
NREFFENS     = TypeVar('Int')
class CONN(CONN_BASE,Generic[MANIFOLD_DIM,NFENS,NREFFENS]):
    
    def __init__(self, manifold_dim: MANIFOLD_DIM, nfens: NFENS, nreffens: NREFFENS):
        self._fens        = None
        self.MANIFOLD_DIM = manifold_dim
        self.NFENS        = nfens
        self.NREFFENS     = nreffens

   
    def local_index(self,fen):
        
        for j in range(0,self.NFENS):
            if ( fen == self._fens(j) ):
                return j
        
        return -1


    def manifold_dim (self):
        return self.MANIFOLD_DIM
    
    def nfens (self):
        return self.NFENS
    
    def fen (self,i):
        return self._fens(i)
    
    def nreffens (self):
        return self.NREFFENS
      
    def conn_code (self):
        return (self.MANIFOLD_DIM*1000000 + self.NREFFENS * 1000 + self.NFENS)

    def get_ref_fen(self, conn, ref_fens, indx):
        pass
    
    def clone (self):
        newconn = CONN(self.MANIFOLD_DIM,self.NFENS,self.NREFFENS)
        newconn._fens = self._fens
        return newconn
    