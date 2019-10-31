lib_LIBRARIES = libmodl.a
libmodl_a_SOURCES = ecell_modl_h20.cc algo_modl.cc ecell_modl.cc  mat_elasticity_iso.cc proto_modl.cc mat_elasticity.cc  ecell_modl_t10.cc ecell_modl_h8.cc 

always_built_SUBDIRS = 

SUBDIRS = $(always_built_SUBDIRS) . 

BOPT = O
PETSC_DEP_LIB = $(PETSC_DEP_LIB)

AM_CXXFLAGS = $(PETSC_DEF) $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated

