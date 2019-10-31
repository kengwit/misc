bindir=./bin/$(ARCH_DIR)
bin_PROGRAMS = a.out 
a_out_SOURCES = main.cc 

EXTRA_DIST = acool.ebc acool.fen acool.fpar acool-gg1.conn acool.hexvue acool.view adapt.sh bigq.ebc bigq.fen bigq.fpar bigq-gg1.conn bigq.hexvue bigq.view cool3.ebc cool3.fen cool3.fpar cool3-gg1.conn cool3.hexvue cool.ebc cool.egf cool.fen cool.fpar cool-gg1.conn cool.hexvue cool.view famuls.ebc famuls.egf famuls.fen famuls.fpar famuls-gg1.conn famuls.hexvue generic.hexvue ghw.sh hexvue.sh junk junk.in junk.info ldomain2.ebc ldomain2.fen ldomain2.fpar ldomain2-gg1.conn ldomain2.hexvue ldomain2.in ldomain.ebc ldomain.fen ldomain.fpar ldomain-gg1.conn ldomain.hexvue ldomain.in ldomain.view lhs main.cc Makefile.am Makefile.in q.ebc q.egf q.fen q.fpar q-gg1.conn q.hexvue q.info q.view rhs runs spiral2.ebc spiral2.fen spiral2.fpar spiral2-gg1.conn spiral2.hexvue spiral.ebc spiral.egf spiral.fen spiral.fpar spiral-gg1.conn spiral.hexvue spiral.view sqcrack.ebc sqcrack.egf sqcrack.fen sqcrack.fpar sqcrack-gg1.conn sqcrack.hexvue sqcrack.hexvue.sh sqcrack.info sqcrack.view test.sh

INCLUDES = @INCLUDES@ -I../src
LDADD = ../../../lib/$(ARCH_DIR)/libsteady_diffusion.a ../../../lib/$(ARCH_DIR)/libfamuls.a 

always_built_SUBDIRS =

SUBDIRS = $(always_built_SUBDIRS) . 

include $(PETSC_MAKE_INCLUDE)
BOPT = g
PETSC_DEP_LIB = $(PETSC_DEP_LIB)

AM_CXXFLAGS = $(PETSC_DEF) $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated

