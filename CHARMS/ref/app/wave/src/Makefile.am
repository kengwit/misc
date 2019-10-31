lib_LIBRARIES = libwave.a
libwave_a_SOURCES = algo_wave.cc ecell_wave.cc ecell_wave_h8.cc proto_wave.cc  ecell_errest_wave_h8.cc ecell_wave_traction_l2.cc ecell_wave_traction_q4.cc ecell_wave_iso_h8.cc ecell_wave_h8_opt.cc ecell_wave_h8_damp.cc ecell_wave_t4_opt.cc ecell_wave_traction_t3.cc ecell_errest_wave_t4.cc mat_viscel.cc mat_viscel_iso.cc

always_built_SUBDIRS = 

SUBDIRS = $(always_built_SUBDIRS)

INCLUDES = @INCLUDES@ -I../../elasticity/src 
LDADD = ../../../lib/$(ARCH_DIR)/libelastcity.a ../../../lib/$(ARCH_DIR)/libfamuls.a ../../../lib/$(ARCH_DIR)/libwave.a

include $(PETSC_MAKE_INCLUDE)
BOPT = O
PETSC_DEP_LIB = $(PETSC_DEP_LIB)
AM_CXXFLAGS = $(PETSC_DEF) $(GIFACE_DEF) -ansi -pedantic -Wno-deprecated

