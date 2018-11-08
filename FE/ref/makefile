LIB_DIR = /cygdrive/c/Users/KWL/Desktop/Research/libraries
BOOST_INC = $(LIB_DIR)/boost_1_56_0
EIGEN_INC = $(LIB_DIR)/eigen-3.2.4
#EIGEN_INC = /home/KWL/libraries/eigen-3.2.1

#LIB_DIR = /Users/kengwit/Desktop/research/libraries
#BOOST_INC = $(LIB_DIR)/boost_1_57_0


CUR_INC = ./
MECHANICS_INC = ./Mechanics
ELEMENT_INC = ./Element
MATERIAL_INC = ./Material

INC = -I$(MECHANICS_INC) -I$(ELEMENT_INC) -I$(MATERIAL_INC) -I$(BOOST_INC) -I$(EIGEN_INC)
GCC = g++

PROGNAME = main

MAIN = main_CPE4F

MYFLAGS = -O3

MECHOBJS = $(MECHANICS_INC)/Globals.o\
		   $(MECHANICS_INC)/Mechanics.o\
		  
ELEMOBJS = $(ELEMENT_INC)/Element.o\
		   $(ELEMENT_INC)/Element_CPE3.o\
		   $(ELEMENT_INC)/Element_CPE4.o\
		   $(ELEMENT_INC)/Element_CPE4F.o\

MATOBJS = $(MATERIAL_INC)/Material.o\
		  $(MATERIAL_INC)/LinearElastic.o\
		  $(MATERIAL_INC)/NeoHookean.o\

MAINOBJ  = $(MAIN).o

OBJS = $(MAINOBJ) $(MECHOBJS) $(ELEMOBJS) $(MATOBJS)

	
$(PROGNAME): $(OBJS)
	$(GCC) $(MYFLAGS) -o $(PROGNAME) $(OBJS)

$(MAIN).o: $(MAIN).cpp 
	$(GCC) $(MYFLAGS) -c $(MAIN).cpp $(INC) 

$(MECHANICS_INC)/%.o: $(MECHANICS_INC)/%.cpp $(MECHANICS_INC)/%.h 
	$(GCC) $(MYFLAGS) -c $< -o $@ $(INC) 
	
$(ELEMENT_INC)/%.o: $(ELEMENT_INC)/%.cpp $(ELEMENT_INC)/%.h
	$(GCC) $(MYFLAGS) -c $< -o $@ $(INC) 
	
$(MATERIAL_INC)/%.o: $(MATERIAL_INC)/%.cpp $(MATERIAL_INC)/%.h
	$(GCC) $(MYFLAGS) -c $< -o $@ $(INC) 

clean:
	rm -f *~ *.o $(OBJS) $(PROGNAME)

