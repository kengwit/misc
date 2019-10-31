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
#ifndef PROTO_FIELD_TRANSFER_H
# define PROTO_FIELD_TRANSFER_H

#include <list>
#include <vector>
#include "algo.h"
#include "gmesh.h"
#include "field.h"
#include "evalpt.h"


template <int NUM_COMPONENTS>
class PROTO_FIELD_TRANSFER {
  
 public:
  
  PROTO_FIELD_TRANSFER () {}
  ~PROTO_FIELD_TRANSFER () {}
  
  /**
     Transfer the src to the dest field.
  */
  void transfer (FIELD<NUM_COMPONENTS> *src, FIELD<NUM_COMPONENTS> *dest)
    {
      FIXED_VECTOR<NUM_COMPONENTS> z(0);
      dest->setto (z); // Zero out 
      // Now start the transfer, from the bottom to the top
      sizet maxlevel = max_level (dest);
      for (sizet level = 0; level < maxlevel+1; level++) {
        const sizet npairs = dest->npairs();
        for (sizet j = 0; j < npairs; j++) {
          FIELD_PAIR<NUM_COMPONENTS> *fp = dest->jth_field_pair (j);
          BFUN *bfun = fp->bfun();
          if (bfun->level() == level) {
            FEN *fen = bfun->fen();
            GCELL *gcell = bfun->gcell(0); // any gcell would do
            // Evaluate basis function values at level lower than the level of bfun
            EVALPT destevalpt (dest, gcell, fen);
            destevalpt.eval ();
            // Evaluate the source field
            FIXED_VECTOR<NUM_COMPONENTS> srcv = src->evaluate (fen, gcell);
            // Loop over all basis functions in the dest field at destevalpt
            double N_J = 1;  // Adjust for PU 
            sizet nbfuns = destevalpt.nbfuns();
            for (sizet J = 0; J < nbfuns; J++) {
              BFUN *dbfun = destevalpt.bfun(J);
              if (dbfun != bfun) {
                double N = destevalpt.N(J);
                BFUN_DOFPARAM_PAIR_ID destdpid = dest->dofparam_id (dbfun->fen());
                CHECK (destdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
                FIXED_VECTOR<NUM_COMPONENTS> dp = dest->field_pair(destdpid)->dofparam();
                srcv.add (-N, dp);
              } else {
                N_J = destevalpt.N(J);
              }
            }
            CHECK(N_J != 0, EXCEPTION_BAD_VALUE,;); // RAW: expensive, use only for debug
            srcv.scale(1/N_J); // Adjust for PU 
            fp->set_dofparam_all (srcv);
          } // if
        } // for
      } // for
      /*
      // if ( dest->npairs() != src->npairs()) cout <<"npairs are not equal"<<"\n";//RAW
      for (sizet i=0; i<dest->npairs(); i++) { //RAW: only for checking purpose
      FIXED_VECTOR<NUM_COMPONENTS> df =  dest->evaluate(dest,dest->bfun(i)->fen()); //RAW 
      BFUN_DOFPARAM_PAIR_ID src_id = src->dofparam_id (dest->bfun(i)->fen());
      if (src_id != INVALID_BFUN_DOFPARAM_PAIR_ID) {
      FIXED_VECTOR<NUM_COMPONENTS> ds =  src->evaluate(src,src->bfun(src_id)->fen());  //RAW
      df.add(-1.0,ds);
      double norm_df = df.l2_norm();
      if (norm_df >DBL_EPSILON)
      cout<<"****error in field transfer:norm(df-ds)=\t" << norm_df<<"\n";//RAW
      }
      }//RAW
      */
    }
/* 
  void transfer (double sval, FIELD<NUM_COMPONENTS> *dest)
    {
      FIXED_VECTOR<NUM_COMPONENTS> z(0);
      dest->setto (z); // Zero out 
      // Now start the transfer, from the bottom to the top
      sizet maxlevel = max_level (dest);
      //dest->debug_display(); //RAW:ashi
      for (sizet level = 0; level < maxlevel+1; level++) {
        const sizet npairs = dest->npairs();
        for (sizet j = 0; j < npairs; j++) {
          FIELD_PAIR<NUM_COMPONENTS> *fp = dest->jth_field_pair (j);
          BFUN *bfun = fp->bfun();
          // if (bfun->fen()->id() == 1708) { //RAW ashi
          //  cout<<"bfun found"<<"\n"; //RAW ashi
         //  } //RAW ashi
          if (bfun->level() == level) {
            FEN *fen = bfun->fen();
            GCELL *gcell = bfun->gcell(0); // any gcell would do
            // Evaluate basis function values at level lower than the level of bfun
            EVALPT destevalpt (dest, gcell, fen);
            destevalpt.eval ();
          ///    if (bfun->fen()->id() == 47) { //RAW ashi 
          ///           destevalpt.debug_display(); 
          ///     } //RAW ashi
            // Evaluate the source field
            FIXED_VECTOR<NUM_COMPONENTS> srcv(sval);
            // Loop over all basis functions in the dest field at destevalpt
            double N_J = 1;  // Adjust for PU 
            sizet nbfuns = destevalpt.nbfuns();
            for (sizet J = 0; J < nbfuns; J++) {
              BFUN *dbfun = destevalpt.bfun(J);
              if (dbfun != bfun) {
                double N = destevalpt.N(J);
                BFUN_DOFPARAM_PAIR_ID destdpid = dest->dofparam_id (dbfun->fen());
                CHECK (destdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
                FIXED_VECTOR<NUM_COMPONENTS> dp = dest->field_pair(destdpid)->dofparam();
                srcv.add (-N, dp);
              } else {
                N_J = destevalpt.N(J);
              }
            }
            CHECK(N_J != 0, EXCEPTION_BAD_VALUE,;); // RAW: expensive, use only for debug
            srcv.scale(1/N_J); // Adjust for PU 
            fp->set_dofparam_all (srcv);
          } // if
        } // for
      } // for
    //   cout << "list of betas*******start**** \n";  
    //   dest->debug_display();
    //   cout << "********end of betas \n";  
    }

*/


void transfer (double sval, FIELD<NUM_COMPONENTS> *dest)
    {
      FIXED_VECTOR<NUM_COMPONENTS> z(0);
      dest->setto (z); // Zero out
      sizet maxlevel = max_level (dest);
      for (sizet level = 0; level < maxlevel+1; level++) {
        const sizet npairs = dest->npairs();
        for (sizet j = 0; j < npairs; j++) {
          FIELD_PAIR<NUM_COMPONENTS> *fp = dest->jth_field_pair (j);
          BFUN *bfun = fp->bfun();
          if (bfun->level() == level) {
            FEN *fen = bfun->fen();
            GCELL *gcell = bfun->gcell(0); // any gcell would do
            EVALPT destevalpt (dest, gcell, fen);
            destevalpt.eval ();
            FIXED_VECTOR<NUM_COMPONENTS> srcv(sval);
            double N_J = 1;  // Adjust for PU
            sizet nbfuns = destevalpt.nbfuns();
            for (sizet J = 0; J < nbfuns; J++) {
              BFUN *dbfun = destevalpt.bfun(J);
              if (dbfun != bfun) {
                double N = destevalpt.N(J);
                BFUN_DOFPARAM_PAIR_ID destdpid = dest->dofparam_id (dbfun->fen());
                CHECK (destdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
                //calculate summation (n*beta*x), it should be equal to x for linear elements and x^2 for quadratic elements 
                FIXED_VECTOR<NUM_COMPONENTS> dp = dest->field_pair(destdpid)->dofparam() * sval;
                srcv.add (-N, dp);
              } else {
                N_J = destevalpt.N(J);
              }
            }
            CHECK(N_J != 0, EXCEPTION_BAD_VALUE,;); // RAW: expensive, use only for debug
            srcv.scale(1/N_J); // Adjust for PU
            fp->set_dofparam_all (srcv);
          } // if
        } // for
      } // for
      // dest->debug_display();
    }
/*
void check_consist  (FIELD<NUM_COMPONENTS> *dest, FIELD<NUM_COMPONENTS> *geometry)
{
      sizet maxlevel = max_level (dest);
      for (sizet level = 0; level < maxlevel+1; level++) {
        const sizet npairs = dest->npairs();
        for (sizet j = 0; j < npairs; j++) {
          FIELD_PAIR<NUM_COMPONENTS> *fp = dest->jth_field_pair (j);
          BFUN *bfun = fp->bfun();
          if (bfun->level() == level) {
            FEN *fen = bfun->fen();
            GCELL *gcell = bfun->gcell(0); // any gcell would do
            EVALPT destevalpt (dest, gcell, fen);
            destevalpt.eval ();
            FIXED_VECTOR<NUM_COMPONENTS> srcv(sval);
            sizet nbfuns = destevalpt.nbfuns();
            for (sizet J = 0; J < nbfuns; J++) {
                double N = destevalpt.N(J);
                BFUN_DOFPARAM_PAIR_ID destdpid = dest->dofparam_id (bfun->fen());
                BFUN_DOFPARAM_PAIR_ID geomdpid = geometry->dofparam_id (bfun->fen());  
                CHECK (destdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
                CHECK (geomdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
                //calculate summation (n*beta*x), it should be equal to x for linear elements and x^2 for quadratic elements 
                FIXED_VECTOR<NUM_COMPONENTS> dp = dest->field_pair(destdpid)->dofparam() * geometry->field_pair(geomdpid)->dofparam();
                srcv.add (N, dp);
            }
            CHECK(N_J != 0, EXCEPTION_BAD_VALUE,;); // RAW: expensive, use only for debug
            srcv.scale(1/N_J); // Adjust for PU
            fp->set_dofparam_all (srcv);
          } // if
        } // for
      } // for
}
 */
/*
   checks if partition of unity holds for all gcells associated with a given fen on n number of random points.
dest = field with PU constants (beta, in this case)
fen_id = PU will be checked at n random points with in gcells supporting fen_id
n = number of random points
*/
void check_pu (FIELD<NUM_COMPONENTS> *dest, int fen_id, int n )
  {
    const sizet npairs = dest->npairs();
    double rint1, rint2;
    int lowest=0, highest=1;
    int range=(highest-lowest)+1; 
    for (sizet j = 0; j < npairs; j++) {
      FIELD_PAIR<NUM_COMPONENTS> *fp = dest->jth_field_pair (j);
      BFUN *bfun = fp->bfun();
      if (bfun->fen()->id() == fen_id) {
        for (sizet i = 0; i < bfun->ngcells(); i++) {
         GCELL *gcell = bfun->gcell(i);
         srand((unsigned)time(0));
         for(int index=0; index< n; index++){
           rint1 = lowest+ (rand()/(RAND_MAX + 1.0)); 
           rint2 = lowest+ (rand()/(RAND_MAX + 1.0));
           if (rint1 < 0) rint1 = -1*rint1;
           if (rint2 < 0) rint2 = -1*rint2;
           POINT param_loc (rint1, rint2); 
           EVALPT evalpt (dest, gcell, param_loc);
           evalpt.eval();
           FIXED_VECTOR<NUM_COMPONENTS> srcv(0);
           for (sizet jj = 0; jj < evalpt.nbfuns(); jj++) {
             BFUN *dbfun = evalpt.bfun(jj);
             double N = evalpt.N(jj);
             BFUN_DOFPARAM_PAIR_ID destdpid = dest->dofparam_id (dbfun->fen());
             CHECK (destdpid != INVALID_BFUN_DOFPARAM_PAIR_ID, EXCEPTION_NULL_PTR,;);
             FIXED_VECTOR<NUM_COMPONENTS> dp = dest->field_pair(destdpid)->dofparam();
             srcv.add (N, dp);
           }
           if ( ((1-srcv(0))) > 1e-6 || ((1-srcv(0)) < -1e-6) ) {
              cout << "****srcv(0) = "<<srcv(0) << "  *****param_loc=  "<<param_loc <<"\n";
         //   CHECK(srcv(0) == 1, EXCEPTION_BAD_VALUE,;);
           }      
         }   
        }
      }
    }
  }
 
 private:
  
  sizet max_level (FIELD<NUM_COMPONENTS> *src) {
    const sizet npairs = src->npairs();
    sizet maxlevel = 0;
    for (sizet j = 0; j < npairs; j++) {
      FIELD_PAIR<NUM_COMPONENTS> *fp = src->jth_field_pair (j);
      sizet level = fp->bfun()->level();
      if (level > maxlevel) maxlevel = level;
    }
    return maxlevel;
  }
  
};

#endif
