from famuls import CHECK
from ref_ctx import REF_CTX
from bfun import BFUN_SET

class ALGO:
    
    def __init__(self, name, mgr ):
        self._name  = name
        self._mgr   = mgr
     
    def name(self):
        return self._name
    
    def mgr(self):
        return self._mgr
    



class ALGO_REFINE(ALGO):
    
    def __init__(self, name, mgr, gmesh ):
        ALGO.__init__(self,name,mgr)
        path = "algorithms/refine/" + self.name()
        #print("path=",end='')
        #print(path)
        
        self._ref_fraction = 1.0
        if (self.mgr().db().param_defined(path + "/ref_fraction")):            
            self._ref_fraction = self.mgr().db().DB_GET_DOUBLE (path + "/ref_fraction")
            
        
        self._h_over_hbar_ref = 1.0
        if (self.mgr().db().param_defined(path + "/h_over_hbar_ref")):            
            self._h_over_hbar_ref = self.mgr().db().DB_GET_DOUBLE(path + "/h_over_hbar_ref");
            
        self._h_over_hbar_unref = 1; 
        if (self.mgr().db().param_defined(path + "/h_over_hbar_unref")):            
            self._h_over_hbar_unref = self.mgr().db().DB_GET_DOUBLE(path + "/h_over_hbar_unref")
        
  
        self._max_ref_level = 100000000 # unlimited
        if (self.mgr().db().param_defined(path + "/max_ref_level")):           
            self._max_ref_level = self.mgr().db().DB_GET_INTEGER (path + "/max_ref_level")
            
        self._gmesh = gmesh
        
        
        self._true_hierarchical = False
        if (self.mgr().db().param_defined (path + "/true_hierarchical")):
            self._true_hierarchical = self.mgr().db().DB_GET_BOOL (path + "/true_hierarchical")

        #print("self._true_hierarchical = ",end='')
        #print(self._true_hierarchical)
            
        self._ref_ctx = REF_CTX(self._gmesh)
   
        self._target_nbfuns = 0 # meaning *not specified*
        if (self.mgr().db().param_defined (path + "/target_nbfuns")):
            self._target_nbfuns = self.mgr().db().DB_GET_INTEGER (path + "/target_nbfuns")

            
    def refinement_strategy(self):
        if (self._true_hierarchical):
            return BFUN_SET.REFINEMENT_STRATEGY.TRUE_HIERARCHICAL
        else:
            return BFUN_SET.REFINEMENT_STRATEGY.QUASI_HIERARCHICAL
        