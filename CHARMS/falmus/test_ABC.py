from abc import ABC, abstractmethod
from point import POINT
from copy import deepcopy

class GCELL(ABC):
    def __init__(self):
        pass    
    
    
    @abstractmethod
    def parent(self):
        pass
    
    def level(self):
        l = 0
        gcell = self        
        while ( gcell.parent() is not None ):
            
            l += 1 
            gcell = gcell.parent()
  
        return l

class GCELL_HEX8(GCELL):
    
    
    NCHILDREN = 8
    
    _child_map  = [ [ [-1,-1,-1], [0,0,0] ],
              [ [ 0,-1,-1], [1,0,0] ],
              [ [ 0, 0,-1], [1,1,0] ],
              [ [-1, 0,-1], [0,1,0] ],
              [ [-1,-1, 0], [0,0,1] ],
              [ [ 0,-1, 0], [1,0,1] ],
              [ [ 0, 0, 0], [1,1,1] ],
              [ [-1, 0, 0], [0,1,1] ] ]


    def __init__(self,parent):
        self._id = -1
        self._parent = parent
        self._child = []
        self._conn = []

    def conn(self):
        return self._conn
        
    def parent(self):
        return self._parent
    
    def map_to_child(self, param_loc, refchild, child_param_loc):
            
            # note: modification of child_param_loc in here, is reflected
            # also outside this function since child_param_loc is an instance
            # and therefore acts like a "reference"
            
            #if (self.nchildren() == 0):
            #    return False
    
            xi    = param_loc(0)
            eta   = param_loc(1)
            theta = param_loc(2)        
            
            if (xi < 0):
                
                if (eta < 0):
                    
                    if (theta < 0):
                        child = self._child[0]
                        child_param_loc[0] = 2 * xi    + 1
                        child_param_loc[1] = 2 * eta   + 1
                        child_param_loc[2] = 2 * theta + 1
                        print('child1')
                    else:
                        child = self._child[4]
                        child_param_loc[0] = 2 * xi    + 1
                        child_param_loc[1] = 2 * eta   + 1
                        child_param_loc[2] = 2 * theta - 1
                        print('child5')
                
                else: # if eta > 0
                    
                    if (theta < 0):
                        child = self._child[3]
                        child_param_loc[0] = 2 * xi    + 1
                        child_param_loc[1] = 2 * eta   - 1
                        child_param_loc[2] = 2 * theta + 1
                        print('child4')

                    else:
                        child = self._child[7]
                        child_param_loc[0] = 2 * xi    + 1
                        child_param_loc[1] = 2 * eta   - 1
                        child_param_loc[2] = 2 * theta - 1
                        print('child8')
            
            else: # if xi > 0
                
                if (eta < 0):
                    
                    if (theta < 0):
                        
                        child = self._child[1]
                        child_param_loc[0] = 2 * xi    - 1
                        child_param_loc[1] = 2 * eta   + 1
                        child_param_loc[2] = 2 * theta + 1
                        print('child2')                    
                        
                    else:
                        
                        child = self._child[5]
                        child_param_loc[0] = 2 * xi    - 1
                        child_param_loc[1] = 2 * eta   + 1
                        child_param_loc[2] = 2 * theta - 1
                        print('child6')                    
          
                else: # if eta > 0
                    
                    if (theta < 0):
                        child = self._child[2]
                        child_param_loc[0] = 2 * xi    - 1
                        child_param_loc[1] = 2 * eta   - 1
                        child_param_loc[2] = 2 * theta + 1
                        print('child3')                    
                        
                    else:
                        
                        refchild[0] = self._child[6]
                        child_param_loc[0] = 2 * xi    - 1
                        child_param_loc[1] = 2 * eta   - 1
                        child_param_loc[2] = 2 * theta - 1
                        #print(refchild[0]._id)
                        #print('child7')                    
      
            
            return True
        
    
    def child(self,j):
        return self._child[j]
        
    def map_to_parent(self):
        
        if ( self._parent is not None ):
            
            for j in range(0,GCELL_HEX8.NCHILDREN):                
                if ( self == self._parent.child(j) ) :
                    return True
                
        
        return False # no parent
            

h0 = GCELL_HEX8(None)
h1 = GCELL_HEX8(h0)
h2 = GCELL_HEX8(h1)
h3 = GCELL_HEX8(h2)

assert(h0.level()==0)
assert(h1.level()==1)
assert(h2.level()==2)
assert(h3.level()==3)
assert(isinstance(h2,GCELL)==True)


param_loc = POINT()
child_param_loc = POINT()
c1 = GCELL_HEX8(h0) 
c2 = GCELL_HEX8(h0)
c3 = GCELL_HEX8(h0)
c4 = GCELL_HEX8(h0)
c5 = GCELL_HEX8(h0)
c6 = GCELL_HEX8(h0)
c7 = GCELL_HEX8(h0)
c8 = GCELL_HEX8(h0)

c1._id = 1
c2._id = 2
c3._id = 3
c4._id = 4
c5._id = 5
c6._id = 6
c7._id = 7
c8._id = 8

c1._parent = h0
c2._parent = h0
c3._parent = h0
c4._parent = h0
c5._parent = h0
c6._parent = h0
c7._parent = h0
c8._parent = h0

h0._child = [c1,c2,c3,c4,c5,c6,c7,c8]

temp=GCELL_HEX8(h0)
ref_temp=[temp]
print('Chec isinstance(ref_temp,list) = ', end='')
print(isinstance(ref_temp,list))
h0.map_to_child(param_loc, ref_temp, child_param_loc)
temp = ref_temp[0]
print('Check (temp._id == c7._id): ',end='')
print(temp._id == c7._id)
print('Check (temp == c7): ',end='')
print(temp == c7)

assert(c1.map_to_parent()==True)
assert(h0.map_to_parent()==False)
