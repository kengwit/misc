from bfun_dofparam_pair_id import INVALID_BFUN_DOFPARAM_PAIR_ID

class DUMMY:
    
    def __init__(self):
        self._bfuns = []
        self._bfun_dofparam_ids = []
        self._nbfuns_total = 0
        self._mem_raise_size = 5
        
        
    def insert_bfun (self,bfun):
        
        # be careful here, it is not simply appending
        if ( len( self._bfuns ) <= self._nbfuns_total ):
            start = len( self._bfuns )
            end   = start + self._mem_raise_size
            for k in range(start, end):
                self._bfuns.append(None)
  
        self._bfuns[ self._nbfuns_total ] = bfun
        self._nbfuns_total += 1
        
        
        if ( len( self._bfun_dofparam_ids ) <= bfun ):
            
            p_size = len(self._bfun_dofparam_ids) # size before resizing
            
            new_total_size = bfun+self._mem_raise_size
            while True:
                self._bfun_dofparam_ids.append(None)
                if (len(self._bfun_dofparam_ids)==new_total_size):
                    break
                        
            for k in range(p_size,len(self._bfun_dofparam_ids)):
                self._bfun_dofparam_ids[k]= INVALID_BFUN_DOFPARAM_PAIR_ID
        
        self._bfun_dofparam_ids[bfun] = self._nbfuns_total-1
   
    
        
            
    def swap(self,temp_loc1, temp_loc2):
        
        #temp_bfun   = self._bfuns[temp_loc1]
        #self._bfuns[temp_loc1] = self._bfuns[temp_loc2]
        #self._bfuns[temp_loc2] = temp_bfun
        
        temp_dpid = self._bfun_dofparam_ids[ temp_loc1 ]
        
        self._bfun_dofparam_ids[  temp_loc1 ] = \
            self._bfun_dofparam_ids[ temp_loc2 ] 
            
        self._bfun_dofparam_ids[ temp_loc2 ] = temp_dpid
          
        print('-----------------inside-----------------')
        print(self._bfun_dofparam_ids)
        print('----------------------------------------')
        
        
        
d = DUMMY()



d.insert_bfun(9)
print(d._bfuns)
print(d._bfun_dofparam_ids)
print( len(d._bfun_dofparam_ids) )

d.swap(0, 9)


#d.insert_bfun(16)
#d.swap(16, d._nbfuns_active)
#print(d._bfuns)
#print(d._bfun_dofparam_ids)
#print( len(d._bfun_dofparam_ids) )
#
#
#d.insert_bfun(18)
##d.swap(18, d._nbfuns_active)
#print(d._bfuns)
#print(d._bfun_dofparam_ids)
#print( len(d._bfun_dofparam_ids) )