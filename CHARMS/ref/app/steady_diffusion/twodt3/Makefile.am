bindir=./bin/$(ARCH_DIR)
bin_PROGRAMS = a.out
a_out_SOURCES = main.cc

EXTRA_DIST = ldomain2.ebc ldomain2.fen ldomain2.fpar ldomain2-gg1.conn ldomain2.hexvue ldomain2.in ldomain.ebc ldomain.fen ldomain.fpar ldomain-gg1.conn ldomain.hexvue ldomain.in ldomain.view main.cc Makefile.am runs test.sh ucsd.ebc ucsd.fen ucsd.fpar ucsd-gg1.conn ucsd.hexvue ucsd.in

INCLUDES = @INCLUDES@ -I../src
LDADD = ../../../lib/$(ARCH_DIR)/libsteady_diffusion.a ../../../lib/$(ARCH_DIR)/libfamuls.a 

always_built_SUBDIRS =

SUBDIRS = $(always_built_SUBDIRS) . 

include $(PETSC_MAKE_INCLUDE)
BOPT = g
PETSC_DEP_LIB = $(PETSC_DEP_LIB)

AM_CXXFLAGS = $(PETSC_DEF)  $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated
