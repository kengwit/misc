
from elements.element_CPE4 import Element_CPE4
from elements.element_C3D10M import Element_C3D10M

def ElementFactory(FE_Type):
    if ( "CPE4" in FE_Type ):
        return Element_CPE4()
    elif ( "C3D10M" in FE_Type ):
        return Element_C3D10M()	
    else:
        raise Exception("Element Factory Fails - Element Type Not Defined\n")
        #print("Element Factory Fails - Element Type Not Defined\n")    
        #return None



 
