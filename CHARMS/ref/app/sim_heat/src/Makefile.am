lib_LIBRARIES = libsim_heat.a
libsim_heat_a_SOURCES = algo_sim_heat.cc


always_built_SUBDIRS = 

SUBDIRS = $(always_built_SUBDIRS) 

INCLUDES = @INCLUDES@ -I../../steady_diffusion/src 
LDADD = ../../../lib/$(ARCH_DIR)/libsteady_diffusion.a ../../../lib/$(ARCH_DIR)/libfamuls.a ../../../lib/$(ARCH_DIR)/libsim_heat.a

include $(PETSC_MAKE_INCLUDE)
BOPT = O
PETSC_DEP_LIB = $(PETSC_DEP_LIB)
AM_CXXFLAGS = $(PETSC_DEF) $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated

