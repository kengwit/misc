bindir=./bin/$(ARCH_DIR)
bin_PROGRAMS = a.out
a_out_SOURCES = main.cc 

EXTRA_DIST = 


INCLUDES = @INCLUDES@ -I..
LDADD = ../../../lib/$(ARCH_DIR)/libelasticity.a ../../../lib/$(ARCH_DIR)/libfamuls.a 
AM_LDFLAGS = 

always_built_SUBDIRS =

SUBDIRS = $(always_built_SUBDIRS) . 

include $(PETSC_MAKE_INCLUDE)
BOPT = O
PETSC_DEP_LIB = $(PETSC_DEP_LIB)

AM_CXXFLAGS = -I../src $(PETSC_DEF) $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated

