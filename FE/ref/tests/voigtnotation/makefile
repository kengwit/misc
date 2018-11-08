LIB_DIR = /Users/kengwit/Desktop/research/libraries
BOOST_INC = $(LIB_DIR)/boost_1_57_0
EIGEN_INC = $(LIB_DIR)/eigen-3.2.4
BLITZ_INC = $(LIB_DIR)/blitz-0.10

LIBS = 
#-lm -lblas -llapack -lopenvdb -ltbb -ltbbmalloc -ltbbmalloc_proxy -lHalf

GCC = g++

PROGNAME = main

MAIN = test

MYFLAGS = -O3

OBJS = $(MAIN).o
			
$(PROGNAME): $(OBJS)
	$(GCC) $(MYFLAGS) -o $(PROGNAME) $(OBJS) $(LIBS)

$(MAIN).o: $(MAIN).cpp 
	$(GCC) $(MYFLAGS) -c $(MAIN).cpp -I$(BOOST_INC) -I$(EIGEN_INC) -I$(BLITZ_INC)


clean:
	rm -f *~ *.o $(PROGNAME)

