# Usage:
# nmake -f makefile 
# nmake -f makefile clean


GCC = g++ 
MYFLAGS = -std=c++11 

PROGNAME = test

$(PROGNAME): *.o
	$(GCC) $(MYFLAGS) -o $(PROGNAME) *.o
	
$(PROGNAME).o: $(PROGNAME).cpp
	$(GCC) -c $(PROGNAME).cpp 
	
clean:
	rm -f *~ *.o $(PROGNAME)
	
	
