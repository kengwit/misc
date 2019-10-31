/*--------------------------------------------------------------------------
--                                                                        --
--         (F)ramework for (A)daptive (MUL)tiphysics (S)imulations        --
--            Copyright (C) 2001, Petr Krysl (pkrysl@ucsd.edu).           --
--                                                                        --
--                 [Pronounce `famulus': scholar's helper]                --
--                                                                        --
--                                                                        --
--  This program is free software; you can redistribute it and/or modify  --
--  it under the terms of the GNU General Public License as published by  --
--  the Free Software Foundation; either version 2 of the License, or     --
--  (at your option) any later version.                                   --
--                                                                        --
--------------------------------------------------------------------------*/
#include "famuls.h"
#include "field.h"
#include "evalpt.h"
#include "square_fixed_matrix.h"
#include "gcell_solid_t4.h"
#include "ecell_wave_t4_opt.h"

const char ECELL_WAVE_T4_OPT::_type_name[] = "wave_t4_opt";


ECELL_WAVE_T4_OPT::ECELL_WAVE_T4_OPT (GCELL *gcell)
  : ECELL_WAVE::ECELL_WAVE (gcell)
{
  _lambda = 0;
  _mu     = 0;
  _rho    = 0;
  { //  for each quadrature point
    _evalpt = 0;
  }
}


bool 
ECELL_WAVE_T4_OPT::assem_internal_forces ()
{
  FIELD_VECTOR *u = _proto->u();
  FIELD_VECTOR *geometry = _proto->geometry();
  
  //this ecell is isotropic and homogeneous so lambda and mu will remain constant.
 
  const double l2mu   = (_lambda + 2*_mu);

  { //  for each quadrature point
    POINT param_loc (0.25, 0.25, 0.25);
    const int nbfuns = _evalpt->nbfuns();
    double sum_N = 0;
    for (sizet i= 0; i<nbfuns; i++) { //RAW
      sum_N = sum_N + _evalpt->N(i); //RAW
    } //RAW
    if (fabs(sum_N)-1 >  1e-7) { //RAW
      cout << "*********tolerance = \t"<< fabs(sum_N)-1<<"\t nbfuns\t"<<nbfuns<<"\n"; //RAW
      //CHECK(fabs(sum_N)-1 < 1e-7,EXCEPTION_BAD_VALUE,;); //RAW  
    }
    // Evaluate derivatives, retrieve displacements, compute strains
    double strain[6]={0,0,0,0,0,0};
    for (int I = 0; I < nbfuns; I++) {
      double N_x = _evalpt->N_x (I);
      double N_y = _evalpt->N_y (I);
      double N_z = _evalpt->N_z (I);
     
      FIXED_VECTOR<3> u_I = u->field_pair (_evalpt->dofparam_id(I))->dofparam ();
      strain[0] += N_x * u_I(0);
      strain[1] += N_y * u_I(1);
      strain[2] += N_z * u_I(2);
      strain[3] += N_y * u_I(0) + N_x * u_I(1);
      strain[4] += N_z * u_I(1) + N_y * u_I(2);
      strain[5] += N_z * u_I(0) + N_x * u_I(2);
    }
    double stress[6]={0,0,0,0,0,0};    
    stress[0] = l2mu * strain[0] + _lambda * (strain[1] + strain[2]);
    stress[1] = l2mu * strain[1] + _lambda * (strain[0] + strain[2]);
    stress[2] = l2mu * strain[2] + _lambda * (strain[1] + strain[0]);
    stress[3] = _mu * strain[3];
    stress[4] = _mu * strain[4];
    stress[5] = _mu * strain[5];
    double fact = (1./6.) * _evalpt->detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_x = _evalpt->N_x (I);
      double N_y = _evalpt->N_y (I);
      double N_z = _evalpt->N_z (I);
      double f_I[3];
      f_I[0] = -fact * (stress[0] * N_x + stress[3] * N_y + stress[5] * N_z);
      f_I[1] = -fact * (stress[1] * N_y + stress[3] * N_x + stress[4] * N_z);
      f_I[2] = -fact * (stress[2] * N_z + stress[4] * N_y + stress[5] * N_x);
      _proto->assemble_f (_evalpt->dofparam_id(I), f_I);
    }
  } //  for each quadrature point
  
  return true;
}

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
bool 
ECELL_WAVE_T4_OPT::assem_k_d_bar ()
{
  FIELD_VECTOR *u = _proto->u();
  FIELD_VECTOR *v = _proto->v();
  FIELD_VECTOR *geometry = _proto->geometry();
  double dt = _proto->dt();
  //this ecell is isotropic and homogeneous so lambda and mu will remain constant.
 
  const double l2mu   = (_lambda + 2*_mu);

  for (int qp = 0; qp < _nqpts; qp++) { //  for each quadrature point
    POINT param_loc (0.25, 0.25, 0.25);
    const int nbfuns = _evalpt->nbfuns();
    // Evaluate derivatives, retrieve displacements, compute strains
    double strain[6]={0,0,0,0,0,0};
    for (int I = 0; I < nbfuns; I++) {
      double N_x = _evalpt->N_x (I);
      double N_y = _evalpt->N_y (I);
      double N_z = _evalpt->N_z (I);
      FIXED_VECTOR<3> u_I = u->field_pair (_evalpt->dofparam_id(I))->dofparam ();
      FIXED_VECTOR<3> v_I = v->field_pair (_evalpt->dofparam_id(I))->dofparam ();
      FIXED_VECTOR<3> acc = _proto->get_acc_t (_evalpt->dofparam_id(I));
      u_I.add(dt,v_I,(dt*dt*0.25),acc);
      strain[0] += N_x * u_I(0);
      strain[1] += N_y * u_I(1);
      strain[2] += N_z * u_I(2);
      strain[3] += N_y * u_I(0) + N_x * u_I(1);
      strain[4] += N_z * u_I(1) + N_y * u_I(2);
      strain[5] += N_z * u_I(0) + N_x * u_I(2);
    }
    double stress[6]={0,0,0,0,0,0};    
    stress[0] = l2mu * strain[0] + _lambda * (strain[1] + strain[2]);
    stress[1] = l2mu * strain[1] + _lambda * (strain[0] + strain[2]);
    stress[2] = l2mu * strain[2] + _lambda * (strain[1] + strain[0]);
    stress[3] = _mu * strain[3];
    stress[4] = _mu * strain[4];
    stress[5] = _mu * strain[5];
    double fact = (1./6.) * _evalpt->detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_x = _evalpt->N_x (I);
      double N_y = _evalpt->N_y (I);
      double N_z = _evalpt->N_z (I);
      double f_I[3];
      f_I[0] = -fact * (stress[0] * N_x + stress[3] * N_y + stress[5] * N_z);
      f_I[1] = -fact * (stress[1] * N_y + stress[3] * N_x + stress[4] * N_z);
      f_I[2] = -fact * (stress[2] * N_z + stress[4] * N_y + stress[5] * N_x);
      _proto->assemble_f_kd (_evalpt->dofparam_id(I), f_I);
    }
  } //  for each quadrature point
  
  return true;
}
#endif

double 
ECELL_WAVE_T4_OPT::strain_energy ()
{
  double total_energy  = 0; 
  double detj = 0;
  {
    double stress[6]  = {0,0,0,0,0,0}; 
    double strain[6] = {0,0,0,0,0,0};
    POINT param_loc (0.25, 0.25, 0.25);
    calc_stress (param_loc, stress, strain);
    total_energy += (stress[0]*strain[0]+ stress[1]*strain[1]+stress[2]*strain[2]
      +stress[3]*strain[3]+stress[4]*strain[4]+stress[5]*strain[5])*_evalpt->detJ ()* (1./6.);
  }
  return (total_energy/2);
}

bool 
ECELL_WAVE_T4_OPT::assem_body_load ()
{
  FIELD_VECTOR *u = _proto->u();
  FIELD_VECTOR *geometry = _proto->geometry();
  LOAD_ON_GCELL_BODY *body_load = this->body_load();
  if (body_load) {
    {
      POINT param_loc (0.25, 0.25, 0.25);
      const int nbfuns = _evalpt->nbfuns ();
      POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
      FIXED_VECTOR<3> force_density = body_load->force_density(at);
      double fact = (1./6.) * _evalpt->detJ ();
      // now loop over all nodes 
      for (int I = 0; I < nbfuns; I++) {
        double N_I = _evalpt->N (I);
        double values[3];
        values[0] = force_density(0)*fact*N_I;
        values[1] = force_density(1)*fact*N_I;
        values[2] = force_density(2)*fact*N_I;
        _proto->assemble_f (_evalpt->dofparam_id(I), values);
      }
    }
  }
  return true;
}



static ECELL_WAVE *
make (GCELL *gcell)
{
  return (new ECELL_WAVE_T4_OPT (gcell));
}

bool
ECELL_WAVE_T4_OPT::register_make_func ()
{
  return ECELL_WAVE::register_make_func (make, string (GCELL_SOLID_T4::TYPE_NAME),
                                         _type_name);
}

double 
ECELL_WAVE_T4_OPT::suggested_time_step () 
{
  FIELD_VECTOR *geometry = _proto->geometry();

#define SUBTR(res, a, b) \
   {res[0] = a[0] - b[0]; res[1] = a[1] - b[1]; res[2] = a[2] - b[2];}
#define COND_MIN(a, b) \
   { if (a != 0) a = min(b, a);else a = b; }
#define LENGTH(d)  sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])  

  double C[6][6];
  double speed_of_sound=0,v_avg=0;
  POINT d; double shortest_side = 0;
  GCELL *gcell        = this->gcell();
  
  for (sizet i = 0; i<gcell->conn()->nfens(); i++) {
    POINT at = gcell->conn()->fen(i)->ref_loc();
    double E     =   _mu*(3*_lambda+2*_mu)/(_lambda+_mu);
    double nu    =   _lambda/(2*(_lambda+_mu));
    v_avg = v_avg + sqrt(E / (1 - 2*nu) / _rho);
  }
  speed_of_sound = v_avg/gcell->conn()->nfens();
    
  SUBTR(d, gcell->conn()->fen(1)->ref_loc(), gcell->conn()->fen(0)->ref_loc());
  COND_MIN(shortest_side, LENGTH(d));
  SUBTR(d, gcell->conn()->fen(2)->ref_loc(), gcell->conn()->fen(0)->ref_loc());
  COND_MIN(shortest_side, LENGTH(d));
  SUBTR(d, gcell->conn()->fen(3)->ref_loc(), gcell->conn()->fen(0)->ref_loc());
  COND_MIN(shortest_side, LENGTH(d));

  SUBTR(d, gcell->conn()->fen(2)->ref_loc(), gcell->conn()->fen(1)->ref_loc());
  COND_MIN(shortest_side, LENGTH(d));
  SUBTR(d, gcell->conn()->fen(3)->ref_loc(), gcell->conn()->fen(1)->ref_loc());
  COND_MIN(shortest_side, LENGTH(d));

  SUBTR(d, gcell->conn()->fen(3)->ref_loc(), gcell->conn()->fen(2)->ref_loc());
  COND_MIN(shortest_side, LENGTH(d));

  return (shortest_side / speed_of_sound);                 
}

bool 
ECELL_WAVE_T4_OPT::assem_lumped_mass ()
{
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();
  {
    int nbfuns = _evalpt->nbfuns();
    double fact = _rho * (1./6.) * _evalpt->detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_I = _evalpt->N (I); 
      _proto->assemble_m (_evalpt->bfun(I)->dofparam_id(), N_I*fact);
    }
  }
  return true;
} 

#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION
bool 
ECELL_WAVE_T4_OPT::assem_consistent_mass ()
{
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();
  {
    POINT param_loc (0.25, 0.25, 0.25);
    EVALPT evalpt (u, gcell(), param_loc);
    evalpt.eval (geometry);
    int nbfuns = evalpt.nbfuns();
    double fact = _rho * gr->ip[qp].Wgt * evalpt.detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_I = evalpt.N (I); 
      for (int J = 0; J < nbfuns; J++) {
        double N_J = evalpt.N (J); 
        _proto->assemble_m (evalpt.bfun(I)->dofparam_id(), evalpt.bfun(J)->dofparam_id(), fact*N_I*N_J);
      }
    }
  }
  return true;
}
#endif 

#if defined(IMPLICIT_SOLVER) && IMPLICIT_SOLVER
bool 
ECELL_WAVE_T4_OPT::assem_bt_k(double K[3][3] )
{
  vector <DOFMAPPER_VECTOR *> dv = dofmappers ();
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();
  const int nbfuns = dv.size ();
  for (int qp = 0; qp < gr->npoints; qp++) {
    POINT param_loc (0.25, 0.25, 0.25);
    EVALPT evalpt (u, gcell (), param_loc);
    evalpt.eval (geometry);
    CHECK (evalpt.nbfuns () == nbfuns, EXCEPTION_BAD_VALUE,;);
    POINT at = geometry->evaluate (geometry, this->gcell(), param_loc);
    double C[6][6];
    double lambda = mat()->var("lambda", at);
    double mu = mat()->var ("mu", at);
    double fact = gr->ip[qp].Wgt * evalpt.detJ ()*_proto->dt()*_proto->dt()*0.25; // beta = 0.25 
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_xI = evalpt.N_x (I);
      double N_yI = evalpt.N_y (I);
      double N_zI = evalpt.N_z (I);
      for (int J = I; J < nbfuns; J++) {
        double N_xJ = evalpt.N_x (J);
        double N_yJ = evalpt.N_y (J);
        double N_zJ = evalpt.N_z (J);
        btcb_3d_iso_symm_C (N_xJ, N_yJ, N_zJ, N_xI, N_yI, N_zI, lambda, mu, fact, K);
      }
    }
  }
  return true;
}

bool 
ECELL_WAVE_T4_OPT::assem_m_bt_k (double m_bt_k[3][3])
{
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();
  {
    POINT param_loc (0.25, 0.25, 0.25);
    EVALPT evalpt (u, gcell(), param_loc);
    evalpt.eval (geometry);
    int nbfuns = evalpt.nbfuns();
    double fact = _rho * gr->ip[qp].Wgt * evalpt.detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_I = evalpt.N (I); 
      for (int J = 0; J < nbfuns; J++) {
        double N_J = evalpt.N (J); 
        _proto->assemble_m_bt_k (evalpt.bfun(I)->dofparam_id(), evalpt.bfun(J)->dofparam_id(), fact*N_I*N_J,m_bt_k);
      }
    }
  }
  return true;
} 

#endif


#if defined(CONSISTENT_MASS_ADAPTION) && CONSISTENT_MASS_ADAPTION 
FIXED_VECTOR<3>
ECELL_WAVE_T4_OPT::momentum()
{
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();
  FIELD_VECTOR *v = _proto->v();
  FIXED_VECTOR<3> p(0); 
  {
    POINT param_loc (0.25, 0.25, 0.25);
    EVALPT evalpt (u, gcell(), param_loc);
    evalpt.eval (geometry);
    int nbfuns = evalpt.nbfuns();
    double fact = _rho * gr->ip[qp].Wgt * evalpt.detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_I = evalpt.N (I); 
      FIXED_VECTOR<3> v_I = v->field_pair (evalpt.dofparam_id(I))->dofparam (); 
      p.add( fact*N_I, v_I);
    }
  }
  return p;
}

double
ECELL_WAVE_T4_OPT::kinetic_energy()
{
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();
  FIELD_VECTOR *v = _proto->v();
  double  p = 0; 
  {
    POINT param_loc (0.25, 0.25, 0.25);
    EVALPT evalpt (u, gcell(), param_loc);
    evalpt.eval (geometry);
    int nbfuns = evalpt.nbfuns();
    double fact = _rho * gr->ip[qp].Wgt * evalpt.detJ ();
    // now loop over all nodes 
    for (int I = 0; I < nbfuns; I++) {
      double N_I = evalpt.N (I); 
      FIXED_VECTOR<3> v_I = v->field_pair (evalpt.dofparam_id(I))->dofparam (); 
       for (int J = 0; J < nbfuns; J++) {
        double N_J = evalpt.N (J); 
        FIXED_VECTOR<3> v_J = v->field_pair (evalpt.dofparam_id(J))->dofparam ();
        p += 0.5*fact*N_I*N_J*v_I.l2_norm()* v_J.l2_norm();
      }
    }
  }
  return p;
}
#endif

void
ECELL_WAVE_T4_OPT::set_watch_point_value (WATCH_POINT *watch_point, POINT &param_loc) 
{
  double stress[6] = {0,0,0,0,0,0};
  double strain[6] = {0,0,0,0,0,0};
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();
  FIELD_VECTOR *v = _proto->v();
  list<string> vars      =  watch_point->var();
  FIXED_VECTOR<3> u_vec = u->evaluate(gcell(), param_loc);
  FIXED_VECTOR<3> v_vec = v->evaluate(gcell(), param_loc);
  POINT at = geometry->evaluate (geometry, this->gcell(), param_loc); //RAW: only for debugging
  calc_stress (param_loc, stress,strain);  
  for (list <string>::iterator i=vars.begin(); i != vars.end(); i++) {
    string  var  = (*i);   
    if (var == "ux") {
      watch_point->set_time_value(var,u_vec(0),_proto->time());
    } else if(var == "uy") {
      watch_point->set_time_value(var,u_vec(1),_proto->time());
    } else if(var == "uz") {
      watch_point->set_time_value(var,u_vec(2),_proto->time());
    } else if(var == "vx") {
      watch_point->set_time_value(var,v_vec(0),_proto->time());
    } else if(var == "vy") {
      watch_point->set_time_value(var,v_vec(1),_proto->time());
    } else if(var == "vz") {
      watch_point->set_time_value(var,v_vec(2),_proto->time());
    } else if(var == "sigma_xx") {
      watch_point->set_time_value(var,stress[0],_proto->time());
    } else if(var == "sigma_yy") {
      watch_point->set_time_value(var,stress[1],_proto->time());
    } else if(var == "sigma_zz") {
      watch_point->set_time_value(var,stress[2],_proto->time());
    } else if(var == "sigma_xy") {
      watch_point->set_time_value(var,stress[3],_proto->time());
    } else if(var == "sigma_yz") { 
      watch_point->set_time_value(var,stress[4],_proto->time());
    } else if(var == "sigma_zx") {
      watch_point->set_time_value(var,stress[5],_proto->time());
    } else if(var == "sigma_yx") {
      watch_point->set_time_value(var,stress[3],_proto->time());
    } else if(var == "sigma_zy") {
      watch_point->set_time_value(var,stress[4],_proto->time());
    } else if(var == "sigma_xz") {
      watch_point->set_time_value(var,stress[5],_proto->time());
    } else if(var == "von_mises") {
      double j2 = ((((stress[0]-stress[1])*(stress[0]-stress[1])) + ((stress[1]-stress[2])*(stress[1]-stress[2]))
                    +((stress[2]-stress[0])*(stress[2]-stress[0])))/6) +stress[3]* stress[3] +stress[4]* stress[4]+stress[5]* stress[5];
      watch_point->set_time_value(var,sqrt(j2),_proto->time());
    } else if(var == "vol_change") {
      double trace_strain = strain[0]+strain[1]+strain[2]; 
      watch_point->set_time_value(var,trace_strain,_proto->time());
    } else if(var == "strain_xx") {
      watch_point->set_time_value(var,strain[0],_proto->time());
    } else if(var == "strain_yy") {
      watch_point->set_time_value(var,strain[1],_proto->time());
    } else if(var == "strain_zz") {
      watch_point->set_time_value(var,strain[2],_proto->time());
    } else if(var == "strain_xy") {
      watch_point->set_time_value(var,strain[3],_proto->time());
    } else if(var == "strain_yz ") {
      watch_point->set_time_value(var,strain[4],_proto->time());
    } else if(var == "strain_zx") {
      watch_point->set_time_value(var,strain[5],_proto->time());
    } else if(var == "strain_energy_density") {
      double w = stress[0]*strain[0]+stress[1]*strain[1]+ stress[2]*strain[2]+2*(stress[3]*strain[3]+stress[4]*strain[4]
                                                                                 +stress[5]*strain[5]);
      watch_point->set_time_value(var,w,_proto->time());
    }
  }
}

FIXED_VECTOR<3> 
ECELL_WAVE_T4_OPT::calc_loc (POINT &param_loc)
{
  FIELD_VECTOR *geometry = _proto->geometry();
  return (geometry->evaluate (geometry, gcell(), param_loc));
}

vector<double> 
ECELL_WAVE_T4_OPT::calc_quants (POINT &param_loc, list<string>vars) 
{
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();
  FIELD_VECTOR *v = _proto->v();
  vector<double> result;
  double stress[6] = {0,0,0,0,0,0};
  double strain[6] = {0,0,0,0,0,0};
  FIXED_VECTOR<3> u_vec = u->evaluate(gcell(),param_loc);
  FIXED_VECTOR<3> v_vec = v->evaluate(gcell(),param_loc);
  calc_stress (param_loc, stress, strain);
  for (list<string>::iterator i = vars.begin(); i != vars.end(); i++) {
    string var = (*i);   
    if (var == "ux") {
      result.push_back(u_vec(0));
    } else if(var == "uy") {
      result.push_back(u_vec(1));
    } else if(var == "uz") {
      result.push_back(u_vec(2));
    } else if(var == "vx") {
      result.push_back(v_vec(0));
    } else if(var == "vy") {
      result.push_back(v_vec(1));
    } else if(var == "vz") {
      result.push_back(v_vec(2));
    } else if(var == "sigma_xx") {
      result.push_back(stress[0]);
    } else if(var == "sigma_xy") {
      result.push_back(stress[3]);
    } else if(var == "sigma_xz") {
      result.push_back(stress[5]);
    } else if(var == "sigma_yx") {
      result.push_back(stress[3]);
    } else if(var == "sigma_yy") {
      result.push_back(stress[1]);
    } else if(var == "sigma_yz") {
      result.push_back(stress[4]);
    } else if(var == "sigma_zx") {
      result.push_back(stress[5]);
    } else if(var == "sigma_zy") {
      result.push_back(stress[4]);
    } else if(var == "sigma_zz") {
      result.push_back(stress[2]);
    } else if(var == "von_mises") {
      double j2 = ((((stress[0]-stress[1])*(stress[0]-stress[1])) + ((stress[1]-stress[2])*(stress[1]-stress[2]))
                    +((stress[2]-stress[0])*(stress[2]-stress[0])))/6) +stress[3]* stress[3] +stress[4]* stress[4]+stress[5]* stress[5];
      result.push_back(sqrt(j2));
    } else if(var == "vol_change") {
      double trace_strain = strain[0] + strain[1] + strain[2];
      result.push_back(trace_strain);
    } else if(var == "strain_xx") {
      result.push_back(strain[0]);
    } else if(var == "strain_yy") {
      result.push_back(strain[1]);
    } else if(var == "strain_zz") {
      result.push_back(strain[2]);
    } else if(var == "strain_xy") {
      result.push_back(strain[3]);
    } else if(var == "strain_yz") {
      result.push_back(strain[4]);
    } else if(var == "strain_zx") {
      result.push_back(strain[5]);
    } else if(var == "strain_energy_density") {
      double w = stress[0]*strain[0]+stress[1]*strain[1]+ stress[2]*strain[2]+2*(stress[3]*strain[3]+stress[4]*strain[4]+stress[5]*strain[5]);
      result.push_back(w);
    } else {
      cout <<"Variable not found\t"<<var<<"\n";
    }
  } 
  return result;
}

bool 
ECELL_WAVE_T4_OPT::calc_stress (POINT &param_loc, double stress[6], double strain[6])
{
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();

  for (int I = 0; I < 6; I++) { strain[I] = 0; }

  EVALPT evalpt (u, gcell (), param_loc);
  evalpt.eval (geometry);
  const int nbfuns = evalpt.nbfuns();
  // Evaluate derivatives, retrieve displacements, compute strains
  for (int I = 0; I < nbfuns; I++) {
    double N_x = evalpt.N_x (I);
    double N_y = evalpt.N_y (I);
    double N_z = evalpt.N_z (I);
    FIXED_VECTOR<3> u_I = u->field_pair (evalpt.dofparam_id(I))->dofparam ();
    strain[0] += N_x * u_I(0);
    strain[1] += N_y * u_I(1);
    strain[2] += N_z * u_I(2);
    strain[3] += N_y * u_I(0) + N_x * u_I(1);
    strain[4] += N_z * u_I(1) + N_y * u_I(2);
    strain[5] += N_z * u_I(0) + N_x * u_I(2);
  }
  const double l2mu = (_lambda + 2*_mu);
  stress[0] = l2mu * strain[0] + _lambda * (strain[1] + strain[2]);
  stress[1] = l2mu * strain[1] + _lambda * (strain[0] + strain[2]);
  stress[2] = l2mu * strain[2] + _lambda * (strain[1] + strain[0]);
  stress[3] = _mu * strain[3];
  stress[4] = _mu * strain[4];
  stress[5] = _mu * strain[5];
  
  return true;
}

bool 
ECELL_WAVE_T4_OPT::calc_stress (POINT &param_loc, double stress[6], double strain[6], double detj)
{
  FIELD_VECTOR *geometry = _proto->geometry();
  FIELD_VECTOR *u = _proto->u();

  EVALPT evalpt (u, gcell (), param_loc);
  evalpt.eval (geometry);
  const int nbfuns = evalpt.nbfuns();
  detj  =  evalpt.detJ ();
  // Evaluate derivatives, retrieve displacements, compute strains
  
  for (int I = 0; I < nbfuns; I++) {
    double N_x = evalpt.N_x (I);
    double N_y = evalpt.N_y (I);
    double N_z = evalpt.N_z (I);
    FIXED_VECTOR<3> u_I = u->field_pair (evalpt.dofparam_id(I))->dofparam ();
    strain[0] += N_x * u_I(0);
    strain[1] += N_y * u_I(1);
    strain[2] += N_z * u_I(2);
    strain[3] += N_y * u_I(0) + N_x * u_I(1);
    strain[4] += N_z * u_I(1) + N_y * u_I(2);
    strain[5] += N_z * u_I(0) + N_x * u_I(2);
  }
  const double l2mu = (_lambda + 2*_mu);
  stress[0] = l2mu * strain[0] + _lambda * (strain[1] + strain[2]);
  stress[1] = l2mu * strain[1] + _lambda * (strain[0] + strain[2]);
  stress[2] = l2mu * strain[2] + _lambda * (strain[1] + strain[0]);
  stress[3] = _mu * strain[3];
  stress[4] = _mu * strain[4];
  stress[5] = _mu * strain[5];
  
  return true;
}

void 
ECELL_WAVE_T4_OPT::attach (PROTO_WAVE *proto, string ggpath)
{
  ECELL_WAVE::attach (proto,ggpath);
  DB *db = proto->db();
  string matname = db->DB_GET_STRING (ggpath + "/material");
  string path = "materials/"+ matname;
  
  if (db->param_defined (path +"/lambda")) {
    _lambda  = db->DB_GET_DOUBLE (path + "/lambda");
  } else {
    if (db->param_defined (path + "/E") &&
        db->param_defined (path + "/nu")) {
      double E  = db->DB_GET_DOUBLE (path + "/E");
      double nu = db->DB_GET_DOUBLE (path + "/nu");  
      _lambda = E * nu / (1 + nu) / (1 - 2*nu);
    } else {
      CHECK (0, EXCEPTION_BAD_VALUE,;);
    }
  }
  
  // this ecell is isotropic and homogeneous and time-independent so lambda and mu may be evaluated
  // once and then will remain constant .

  if (db->param_defined (path + "/mu")) {
    _mu = db->DB_GET_DOUBLE (path + "/mu");
  } else {
    if (db->param_defined (path + "/E") &&
        db->param_defined (path + "/nu")) {
      double E  = db->DB_GET_DOUBLE (path + "/E");
      double nu = db->DB_GET_DOUBLE (path + "/nu"); 
      _mu = E / 2 / (1 + nu);
    } else {
      CHECK (0, EXCEPTION_BAD_VALUE,;);
    }
  }
  
  if (db->param_defined (path + "/rho")) {
    _rho = db->DB_GET_DOUBLE (path + "/rho");
  }

  { //  for each quadrature point
    POINT param_loc (0.25, 0.25, 0.25);
    _evalpt = new EVALPT (proto->u(), gcell (), param_loc);
    _evalpt->eval (proto->geometry());
    CHECK(!isnan(_evalpt->detJ ()), EXCEPTION_BAD_VALUE,;);
    for (int I = 0; I < _evalpt->nbfuns(); I++) {
      CHECK(!isnan(_evalpt->N (I)), EXCEPTION_BAD_VALUE,;);
    }
  }
}

ECELL_WAVE_T4_OPT::~ECELL_WAVE_T4_OPT() {
  { //  for each quadrature point
    delete _evalpt;
  }
}
