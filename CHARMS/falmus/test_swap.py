
class DUMMY:
    
    def __init__(self,bfuns):
        
        self._bfuns = bfuns


    def swap(self,temp_loc1, temp_loc2):
        
        temp_bfun   = self._bfuns[temp_loc1]
        self._bfuns[temp_loc1] = self._bfuns[temp_loc2]
        self._bfuns[temp_loc2] = temp_bfun
        
        

bfuns = [11,22,33]

d = DUMMY(bfuns.copy())  
print(d._bfuns)

temp_loc1 = 0
temp_loc2 = 1

d.swap(temp_loc1,temp_loc2)      
assert(d._bfuns[temp_loc1] == bfuns[temp_loc2])
assert(d._bfuns[temp_loc2] == bfuns[temp_loc1])

