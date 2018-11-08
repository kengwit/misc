LIBDIR = /cygdrive/c/cygwin/home/KWL/libraries
EIGENINC = -I$(LIBDIR)/eigen-3.2.1

G++ = g++

PROGNAME = main

MAIN = main

MYFLAGS = -O3

OBJS = $(MAIN).o\
		ParDerGen.o
		  
$(PROGNAME): $(OBJS)
	$(G++) $(MYFLAGS) -o $(PROGNAME) $(OBJS) $(EIGENINC)

$(MAIN).o: $(MAIN).cpp *.h
	$(G++) $(MYFLAGS) -c $(MAIN).cpp $(EIGENINC)

ParDerGen.o: ParDerGen.cpp *.h
	$(G++) $(MYFLAGS) -c ParDerGen.cpp $(EIGENINC)
	

clean:
	rm -f *~ *.o $(PROGNAME)

