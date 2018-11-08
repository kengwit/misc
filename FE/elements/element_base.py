from abc import ABCMeta, abstractmethod


class Element_Base(metaclass=ABCMeta):
    
    def __init__(self,ndim,nnodes,ndof):
        self.id      = 0
        self.nsdm    = 0
        self.ndim    = ndim
        self.nnodes  = nnodes
        self.ndof    = ndof
        
        self.conn    = []
        self.Ulocal  = []
        self.dUlocal = []
        
        self.xref     = None
        self.xcurrent = None
        
        self.nvoigt = 0 
        self.Bmatl = None
        self.Gmatl = None
        self.Amatl = None
        self.Svecl = None
        
        self.Fl    = None
        self.Fl_t  = None
        
        self.ke    = None
        self.fintl = None
        
        self.weights = None
        
    @abstractmethod
    def compute_shape(self):
        pass
    
    def compute_stress(self):
        pass
    
    def assemble(self):
        pass
    
    



    
        
    