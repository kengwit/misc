LIB_DIR = /cygdrive/c/Users/KWL/Desktop/Research/libraries
BOOST_INC = $(LIB_DIR)/boost_1_57_0
EIGEN_INC = $(LIB_DIR)/eigen-3.2.4

LIBS = 
#-lm -lblas -llapack -lopenvdb -ltbb -ltbbmalloc -ltbbmalloc_proxy -lHalf

GCC = g++

PROGNAME = main

MAIN = tensormult

MYFLAGS = -O3 -std=c++11

OBJS = $(MAIN).o
			
$(PROGNAME): $(OBJS)
	$(GCC) $(MYFLAGS) -o $(PROGNAME) $(OBJS) $(LIBS)

$(MAIN).o: $(MAIN).cpp 
	$(GCC) $(MYFLAGS) -c $(MAIN).cpp -I$(BOOST_INC) -I$(EIGEN_INC)


clean:
	rm -f *~ *.o $(PROGNAME)

