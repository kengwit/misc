from point import POINT

v = POINT()

assert(v.size()==3)

v[0]=1
print(v(0))


def func(v):
    # modification possible because v is a list
    v[1]=22.
    
    
func(v)    
print(v(1))