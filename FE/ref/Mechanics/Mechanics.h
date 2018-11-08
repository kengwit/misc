#ifndef MECHANICS_H
#define MECHANICS_H

#include "Globals.h"

#define FORCE_BC 0
#define DISPL_BC 1

void mech_loop_through_elements(int isw);
void mech_implicit_solve();
void mech_solve_linearized_system(double scale_disp_bc);
	

#endif