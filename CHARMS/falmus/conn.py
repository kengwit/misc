from famuls import CHECK
from abc import ABC, abstractmethod
from enum import Enum
from typing import TypeVar, Generic, overload

class CONN_MANIFOLD_DIM(Enum):    
    CONN_0_MANIFOLD = 0
    CONN_1_MANIFOLD = 1
    CONN_2_MANIFOLD = 2
    CONN_3_MANIFOLD = 3
        
   
#T  = TypeVar('CONN_MANIFOLD_DIM')
#NF = TypeVar('Int')
#NR = TypeVar('Int')
#class CONN(Generic[T]):
class CONN():
     
    #def __init__(self,MANIFOLD_DIM: T, NFENS: NR , NREFFENS: NF):
    def __init__(self,MANIFOLD_DIM, NFENS, NREFFENS):
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

    def clone (self):
        newconn = CONN(self.MANIFOLD_DIM,self.NFENS,self.NREFFENS)
        newconn._fens = self._fens
        return newconn

    def get_ref_fen(self, conn, ref_fens, indx):
        
        if ( self.MANIFOLD_DIM == CONN_MANIFOLD_DIM.CONN_0_MANIFOLD and
             self.NFENS == 1 and self.NREFFENS == 1 ):
            
            #print('get_ref_fen(CONN_MANIFOLD_DIM.CONN_0_MANIFOLD,1,1)')
            CHECK (len(ref_fens) == 1)
            CHECK (conn.manifold_dim() == self.manifold_dim())
            CHECK (conn.nfens() == self.nfens())
            CHECK (indx == 0)
            
            if (conn.fen(0) == self.fen(0)):
                return ref_fens[0]
            else:
                return 0
        
        elif ( self.MANIFOLD_DIM == CONN_MANIFOLD_DIM.CONN_1_MANIFOLD and
               self.NFENS == 2 and self.NREFFENS == 1 ):           
        
            #print('get_ref_fen(CONN_MANIFOLD_DIM.CONN_1_MANIFOLD,2,1)')
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
        
        elif ( self.MANIFOLD_DIM == CONN_MANIFOLD_DIM.CONN_2_MANIFOLD and
               self.NFENS == 4 and self.NREFFENS == 1 ):           
            
            CHECK (len(ref_fens) == 1)
            CHECK (conn.manifold_dim() == self.manifold_dim())
            CHECK (conn.nfens() == self.nfens())
            CHECK (indx == 0)
            
            def CHK4(n0,n1,n2,n3):
                return ((conn.fen(0) == self.fen(n0)) and
                        (conn.fen(1) == self.fen(n1)) and
                        (conn.fen(2) == self.fen(n2)) and
                        (conn.fen(3) == self.fen(n3)) )
            
            if CHK4(0,1,2,3):
                return ref_fens[0]
            elif CHK4(1,2,3,0):
                return ref_fens[0]
            elif CHK4(2,3,0,1):
                return ref_fens[0]
            elif CHK4(3,0,1,2):
                return ref_fens[0]
            elif CHK4(3,2,1,0): # and then try the mirror configurations
                return ref_fens[0]
            elif CHK4(2,1,0,3):
                return ref_fens[0]
            elif CHK4(1,0,3,2):
                return ref_fens[0]
            elif CHK4(0,3,2,1):
                return ref_fens[0]
            else:
                return 0 

        elif ( self.MANIFOLD_DIM == CONN_MANIFOLD_DIM.CONN_3_MANIFOLD and
               self.NFENS == 8 and self.NREFFENS == 1 ):           
            
            CHECK (len(ref_fens) == 1)
            CHECK (conn.manifold_dim() == self.manifold_dim())
            CHECK (conn.nfens() == self.nfens())
            CHECK (indx == 0)
            
            # We could proceed as with the quadrilateral, but another option
            # is just to count the number of present nodes
            num_matched = 0
            for j in range(0,self.nfens()):
               to_match = self.fen(j)
               for k in range(0,self.nfens()):             
                   if (to_match == conn.fen(k)):
                       num_matched += 1
                       break
      
            if (num_matched == 8):
                return ref_fens[0]
            else:
                return 0
                        
        else:
            raise NotImplementedError
            
        
        
def CONN_POINT_1():
    return CONN(CONN_MANIFOLD_DIM.CONN_0_MANIFOLD, NFENS = 1, NREFFENS = 1)

def CONN_LINE_2():
    return CONN(CONN_MANIFOLD_DIM.CONN_1_MANIFOLD, NFENS = 2, NREFFENS = 1)

def CONN_SURF_4():
    return CONN(CONN_MANIFOLD_DIM.CONN_2_MANIFOLD, NFENS = 4, NREFFENS = 1)

def CONN_SOLID_8():
    return CONN(CONN_MANIFOLD_DIM.CONN_3_MANIFOLD, NFENS = 8, NREFFENS = 1)
 
    
class CONN_REF:
    
    def __init__(self,conn):
        self._conn     = conn.clone()
        self._ref_fens = []
        
    def get_ref_fen(self, conn, indx):
        return self._conn.get_ref_fen(conn, self._ref_fens, indx);

    # Refine the conn_ref by generating refinement nodes.
    def refine(self,fen_maker):
        
        ####in C++, fen_maker is a pointer, do we want a reference wrapper in Python?
        raise NotImplementedError
        
        
        self._ref_fens = []
        nreffens = self._conn.nreffens();
        if (nreffens > 0):
            # Generate refinement nodes.
            for j in range(0,nreffens):
                fen = fen_maker(j)
                self._ref_fens.append(fen)
        