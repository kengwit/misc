
#from abc import ABC, abstractmethod
#
#class ENUMERATOR(ABC):
#    
#    def _init_(self,container):
#        self._container = container
#
#    @abstractmethod
#    def __iter__(self):
#        pass    
    

class BFUN_SET:
    
    class bfun_enumerator_t():        
        
        def __init__(self,bfuns):
            self._bfuns = bfuns
            
        def __iter__(self):
            for bfun in self._bfuns:
                if ( bfun != -1 ):
                    yield bfun 
    
    class bfun_enumerator_all_t():        
        
        def __init__(self,bfuns):
            self._bfuns = bfuns
            
        def __iter__(self):
            for bfun in self._bfuns:
                yield bfun 
                
                
    def __init__(self,bfuns):
        self._bfuns = bfuns
        
        
    def bfun_enumerator(self):
        e = self.bfun_enumerator_t(self._bfuns)
        return e
    
    def bfun_enumerator_all(self):
        e = self.bfun_enumerator_all_t(self._bfuns)
        return e
        
        

bfuns = [1,-1,10,-1,-2]
t = BFUN_SET(bfuns)
e = t.bfun_enumerator()
for bfun in e:
    print(bfun)

        

            
                