lib_LIBRARIES = libelasticity.a
libelasticity_a_SOURCES = algo_elasticity.cc ecell_elasticity.cc ecell_elasticity_t4.cc ecell_elasticity_h8.cc mat_elasticity_iso.cc proto_elasticity.cc mat_elasticity.cc ecell_errest_elasticity_t4.cc ecell_errest_elasticity_h8.cc load_on_gcell.cc load_on_gcell_body.cc load_on_gcell_traction.cc ecell_elasticity_traction_l2.cc ecell_elasticity_plane_strain_q4.cc ecell_errest_elasticity_plane_strain_q4.cc ecell_elasticity_h8_bbar.cc ecell_elasticity_traction_q4.cc 

always_built_SUBDIRS = 

SUBDIRS = $(always_built_SUBDIRS) . 

include $(PETSC_MAKE_INCLUDE)
BOPT = O
PETSC_DEP_LIB = $(PETSC_DEP_LIB)

AM_CXXFLAGS = $(PETSC_DEF) $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated

