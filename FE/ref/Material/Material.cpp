#include "Material.h"
#include "LinearElastic.h"
#include "NeoHookean.h"

Material * Material::Factory(int choice)
{
	switch (choice) {
	case LINEAR_ELASTIC:
		return new LinearElastic();
		break;
	
	case NEOHOOKEAN:
		return new NeoHookean();
		break;
		
	default:
		cout << "Undefined material type";
		return NULL;
	}

}
	
void MaterialDriver(unsigned int mat, unsigned int el, unsigned int sec_type, unsigned int gp, ViewD2 & eps_t, ViewD2 & eps_tau, ViewD2 & Dmat, ViewD1 & Svec, ViewD1 & svars, bool tangent_flag)
{
	fem::material_database[mat]->constitutive_update(el, gp, sec_type, fem::material_database[mat]->props, eps_t, eps_tau, Dmat, Svec, svars, tangent_flag);	
}

