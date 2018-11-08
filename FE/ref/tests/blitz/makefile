LIB_DIR = /cygdrive/c/Users/KWL/Desktop/Research/libraries
BLITZ_INC = $(LIB_DIR)/blitz-0.10

LIBS = 

GCC = g++

PROGNAME = main

MAIN = tensor

MYFLAGS = -O3 -std=c++11

OBJS = $(MAIN).o
			
$(PROGNAME): $(OBJS)
	$(GCC) $(MYFLAGS) -o $(PROGNAME) $(OBJS) $(LIBS)

$(MAIN).o: $(MAIN).cpp 
	$(GCC) $(MYFLAGS) -c $(MAIN).cpp -I$(BLITZ_INC)

clean:
	rm -f *~ *.o $(PROGNAME)

