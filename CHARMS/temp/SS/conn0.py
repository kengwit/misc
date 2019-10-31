from famuls import CHECK
from abc import ABC, abstractmethod
from enum import Enum

class CONN_MANIFOLD_DIM(Enum):
  CONN_0_MANIFOLD = 0
  CONN_1_MANIFOLD = 1
  CONN_2_MANIFOLD = 2
  CONN_3_MANIFOLD = 3
        


'''
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
        
'''        

class CONN(ABC):
    
    def __init__(self,MANIFOLD_DIM,NFENS,NREFFENS):
        self._fens        = None
        self.MANIFOLD_DIM = MANIFOLD_DIM
        self.NFENS        = NFENS
        self.NREFFENS     = NREFFENS
   
    def local_index(self,fen):
        
        for j in range(0,self.NFENS):
            if ( fen == self._fens[j] ):
                return j
        
        return -1

    def manifold_dim (self):
        return self.MANIFOLD_DIM

    def fen(self,i):
        return self._fens[i]
    
    def nfens(self):
        return self.NFENS

    def nreffens(self):
        return self.NREFFENS
      
    def conn_code (self):
        return (self.MANIFOLD_DIM*1000000 + self.NREFFENS * 1000 + self.NFENS)

    @abstractmethod
    def clone (self):
        pass
        #newconn = CONN(self.MANIFOLD_DIM,self.NFENS,self.NREFFENS)
        #newconn._fens = self._fens
        #return newconn

    @abstractmethod
    def get_ref_fen(self, conn, ref_fens, indx):
        pass
    
    
#    def swap(n0,n1):
#        afen    = fen(n0)
#        fen(n0) = fen(n1)
#        fen(n1) = afen

# ---------------------------------------------------------    
# Connectivity of an isolated node.
#   The refinement is just another isolated node.
# ---------------------------------------------------------    
class CONN_POINT_1(CONN):
    def __init__(self):
        CONN.__init__(self,MANIFOLD_DIM = CONN_MANIFOLD_DIM.CONN_0_MANIFOLD,
                           NFENS = 1,
                           NREFFENS = 1)
    
    def get_ref_fen(self, conn, ref_fens, indx):
        CHECK (len(ref_fens) == 1)
        CHECK (conn.manifold_dim() == self.manifold_dim())
        CHECK (conn.nfens() == self.nfens())
        CHECK (indx == 0)
        
        if (conn.fen(0) == self.fen(0)):
            return ref_fens[0]
        else:
            return 0
        
    def clone(self):
        newconn = CONN_POINT_1()
        newconn._fens = self._fens
        return newconn
        
# ---------------------------------------------------------    
#   Connectivity of a 2-node line.  
#   Numbering of nodes: 
#   0-------------------1
#   The refinement consists of one interior node, and two vertex nodes.
#   X---------0---------X
# ---------------------------------------------------------    
class CONN_LINE_2(CONN):
    def __init__(self):
        CONN.__init__(self,MANIFOLD_DIM = CONN_MANIFOLD_DIM.CONN_1_MANIFOLD,
                           NFENS = 2,
                           NREFFENS = 1)
    
    def get_ref_fen(self, conn, ref_fens, indx):
        
        CHECK (len(ref_fens) == 1)
        CHECK (conn.manifold_dim() == self.manifold_dim())
        CHECK (conn.nfens() == self.nfens())
        CHECK (indx == 0)
        
        if ( ( conn.fen(0) == self.fen(0) ) and 
             ( conn.fen(1) == self.fen(1) ) ):
            return ref_fens[0]
        elif ( ( conn.fen(0) == self.fen(1) ) and 
               ( conn.fen(1) == self.fen(0) ) ):
            return ref_fens[0]
        else:
            return 0
        
    def clone(self):
        newconn = CONN_LINE_2()
        newconn._fens = self._fens
        return newconn        