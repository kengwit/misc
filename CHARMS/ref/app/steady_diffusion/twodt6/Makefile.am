bindir=./bin/$(ARCH_DIR)
bin_PROGRAMS = a.out
a_out_SOURCES = main.cc

EXTRA_DIST = Makefile.am

INCLUDES = @INCLUDES@ -I../src
LDADD = ../../../lib/$(ARCH_DIR)/libsteady_diffusion.a ../../../lib/$(ARCH_DIR)/libfamuls.a 

always_built_SUBDIRS =

SUBDIRS = $(always_built_SUBDIRS) . 

include $(PETSC_MAKE_INCLUDE)
BOPT = g
PETSC_DEP_LIB = $(PETSC_DEP_LIB)

AM_CXXFLAGS = $(PETSC_DEF) $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated