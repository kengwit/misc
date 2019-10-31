bindir=./bin/$(ARCH_DIR)
bin_PROGRAMS = a.out
a_out_SOURCES = main.cc

EXTRA_DIST = ldomain2.ebc ldomain2.fen ldomain2.fpar ldomain2-gg1.conn ldomain2.in ldomain.ebc ldomain.egf ldomain.fen ldomain.fpar ldomain-gg1.conn ldomain.in Makefile.am 


INCLUDES = @INCLUDES@ -I../src
LDADD = ../../../lib/$(ARCH_DIR)/libsteady_diffusion.a ../../../lib/$(ARCH_DIR)/libfamuls.a 

always_built_SUBDIRS =

SUBDIRS = $(always_built_SUBDIRS) . 

include $(PETSC_MAKE_INCLUDE)
BOPT = O
PETSC_DEP_LIB = $(PETSC_DEP_LIB)

AM_CXXFLAGS = $(PETSC_DEF) $(GIFACE_DEF) $(SILO_DEF)  -ansi -pedantic -Wno-deprecated

