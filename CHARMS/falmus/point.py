import numpy as np

class FIXED_VECTOR:

    def __init__(self,n):        
        self._data = np.zeros(n)

    #def __getitem__(self,i):
    #    return self._data[i]

    def __setitem__(self,i,val):
        self._data[i]=val

    def __call__(self,i): 
        return self._data[i]

    def size(self):
        return self._data.size



def POINT():
    return FIXED_VECTOR(3)


def POINT4():
    return FIXED_VECTOR(4)
