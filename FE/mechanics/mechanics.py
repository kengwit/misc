import numpy as np
from enum import Enum

class LOOP(Enum):
    DO_SHAPE    = 0
    DO_STRAIN   = 1
    DO_KMAT_RES = 2
    
class FORMULATION(Enum):
    SMALL_STRAIN       = 0
    UPDATED_LAGRANGIAN = 1
    TOTAL_LAGRANGIAN   = 2 
    
class Mechanics():
        
    def __init__(self):        
        self.ElementMap = {}
        self.X_coords = {}
        self.x_coords = {}
        
        self.nelems  = 0
        self.nnodes  = 0
        self.ndofs   = 0
        self.nqps    = 0
        self.ndim    = 0
        self.matdb   = {}
		
        self.U_tau   = {}
        self.U_t     = {}
        self.dU      = {}
        self.fext    = {}
        self.fint    = {}
        self.frct    = {}
        self.bcdofs  = []
        self.bcvals  = []
        #self.isw     = LoopSwitch(Enum)
        #self.MapNodeIDtoIndex = {}

	  
    def mech_implicit_solve(self):
        for i in range(0,2):
            print('step: %d'%i)
            self.loop_through_elements(LOOP.DO_SHAPE)
            #self.loop_through_elements(LOOP.DO_STRAIN)
            #self.loop_through_elements(LOOP.DO_KMAT_RES)
       

    def loop_through_elements(self,ISW):        
        
        for ie,el in self.ElementMap.items():
            
            if ( ISW is LOOP.DO_SHAPE ):
                
                #print("elem: %d, conn = "%(el.id),end=" ")
                #print(el.conn,end="\n")
                #el.Ulocal  = []
                #el.dUlocal = []
                for i in range(0,len(el.conn)):
                    node_id    = el.conn[i]                    
                    #el.Ulocal.append( self.U_tau[ node_id ] )
                    #el.dUlocal.append( self.dU[ node_id ] )
                    
                    for k in range(0,el.ndim):
                        el.xref[k,i] = self.X_coords[ node_id ][k]
                        
                        
                        
                
                    #print("i=%d, xlocal = "%i)
                    #print(el.xref,end="\n")
                    
                
                el.compute_shape()
                
            if ( ISW is LOOP.DO_STRAIN ):
                pass
                
            if ( ISW is LOOP.DO_KMAT_RES ):
                
                #print("elem: %d, conn = "%(el.id),end=" ")
                #print(el.conn,end="\n")
                
                for i in range(0,len(el.conn)):
                    pass
                    # to do copy history
                
                #print(el.xlocal)
            