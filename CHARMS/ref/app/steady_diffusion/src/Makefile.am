lib_LIBRARIES = libsteady_diffusion.a
libsteady_diffusion_a_SOURCES = algo_steady_diffusion.cc ecell_errest_steady_diffusion_h8.cc ecell_errest_steady_diffusion_l2.cc ecell_errest_steady_diffusion_q4.cc ecell_errest_steady_diffusion_t3.cc ecell_errest_steady_diffusion_t4.cc ecell_errest_steady_diffusion_t10.cc ecell_errest_steady_diffusion_t6.cc ecell_steady_diffusion.cc ecell_steady_diffusion_h8.cc ecell_steady_diffusion_l2.cc ecell_steady_diffusion_q4.cc ecell_steady_diffusion_t3.cc ecell_steady_diffusion_t4.cc ecell_steady_diffusion_t10.cc ecell_steady_diffusion_t6.cc proto_steady_diffusion.cc mat_diffusion_iso.cc mat_diffusion.cc mat_diffusion_aniso.cc ecell_steady_diffusion_t3_aniso.cc 

always_built_SUBDIRS = 

SUBDIRS = $(always_built_SUBDIRS) . 

include $(PETSC_MAKE_INCLUDE)
BOPT = O
PETSC_DEP_LIB = $(PETSC_DEP_LIB)
AM_CXXFLAGS = $(PETSC_DEF) $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated

