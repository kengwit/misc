bindir=./bin/$(ARCH_DIR)
bin_PROGRAMS = a.out
a_out_SOURCES = main.cc 

EXTRA_DIST = 


INCLUDES = @INCLUDES@ -I.. -I../src
LDADD = ../../../lib/$(ARCH_DIR)/libmodl.a ../../../lib/$(ARCH_DIR)/libfamuls.a 
AM_LDFLAGS = 

always_built_SUBDIRS =

SUBDIRS = $(always_built_SUBDIRS) . 

AM_CXXFLAGS = -I../src $(PETSC_DEF) $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated

