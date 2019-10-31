#include "subspaceit.h"
#include <stdio.h>
#include <math.h>
#include "flotmtrx.h"
#include "skyline.h"
#include "flotarry.h"
#include "gjacobi.h"
#include "ooofs_error.h"

#include "logger_stream.h"
#include "logger_stream_mgr.h"


#if !defined(max)
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#if !defined(min)
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

//#define DETAILED_REPORT 1

namespace OOOFS {

SubspaceIteration ::  SubspaceIteration ()  {
  nc     = 0 ;
  n = 0 ;
  nitem = 40;  // max number of iterations
  solved = 0 ;
}


SubspaceIteration ::  ~SubspaceIteration () {
}

bool 
SubspaceIteration::solve (LOGGER_STREAM &ls, SparseMtrx* a, SparseMtrx* b, FloatArray* _eigv, 
                          FloatMatrix* _r, double rtol,
			  int max_iter, int nroot)
{
  FloatArray temp, w, d, tt, rtolv, eigv ;
  FloatMatrix r;
  int nn,nd,nc1,i,j,k,l,ij=0,nite,iconv,is ;
  double rt,art,brt,eigvt,dif;
  FloatMatrix ar, br, vec;
  bool converged = false;

  nitem = max_iter;

  GJacobi *mtd = new GJacobi();
  nc = min(2*nroot,nroot+8) ;
  ls << " dimension of the reduced problem=" << nc << endl_notime;
  //
  // check matrix size
  // 
  if ((!a) || (!b)) _error ("SubspaceIteration :: solveYourselfAt : matrices are not defined\n");
  if (a->giveNumberOfColumns() != b->giveNumberOfColumns())
    _error("SubspaceIteration :: solveYourselfAt : matrices size mismatch\n");
  // check matrix for factorization support
  if (!a->canBeFactorized()) _error ("SubspaceIteration :: a matrix not support factorization");
  //
  // check for wery small problem
  //
  nn = a->giveNumberOfColumns() ;
  if (nc > nn ) nc = nn;
  
  ar.resize (nc,nc); ar.zero();
  br.resize (nc,nc); br.zero();
  
  //
  // creation of initial iteration vectors
  //
  iconv = 0 ;
  nd = nn/nc ;
  nc1 = nc - 1 ;
  /*
    w = new FloatArray (nn);
    d = new FloatArray (nc);
    tt = new FloatArray (nn);
    rtolv= new FloatArray (nc);
  */
  w.resize (nn); w.zero();
  d.resize (nc); d.zero();
  tt.resize (nn); tt.zero();
  rtolv.resize (nc); rtolv.zero();
  vec.resize(nc,nc); vec.zero();  // eigen vectors of reduced problem
  // check matrix for storing resulted eigen vectors at the end
  if (_r == NULL) _error ("solveYourselfAt: unknown eigen vectors mtrx");
  if ((_r->giveNumberOfRows () != nn) || (_r->giveNumberOfColumns() != nroot))
    _error ("solveYourselfAt: _r size mismatch");
  // check array for storing eigenvalues
  if (_eigv == NULL) _error("solveYourselfAt: unknown eigenvalue array");
  if (_eigv->giveSize() != nroot) _error ("solveYourselfAt: eigv size mismatch");
  // 
  // create work arrays
  // 
  /*
    r     = new FloatMatrix(nn,nc);
    eigv  = new FloatArray (nc);
  */
  r.resize (nn, nc); r.zero();
  eigv.resize (nc); eigv.zero();
  
  for (i = 1; i<= nn ; i++) {
    r.at(i,1) = b->at(i,i);
    w.at(i) = b->at(i,i)/a->at(i,i);
  } // label 10
  l = nn - nd;
  for (j = 2; j<= nc ; j++) {
    rt = 0.0 ;
    for (i = 1; i<= l ; i++) 
      if (fabs(w.at(i)) >= rt) { rt = fabs(w.at(i)); ij = i;}
    for (i = l; i<= nn; i++) 
      if (fabs(w.at(i)) > rt) { rt = fabs(w.at(i)); ij = i;}
    tt.at(j) = ij ;
    w.at(ij) = 0. ;
    l -= nd ;
    r.at(ij,j) = 1.0;
  }
# ifdef DETAILED_REPORT
  printf("SubspaceIteration :: solveYourselfAt: Degrees of freedom invoked by initial vectors :\n");
  tt.printYourself(); 
  printf ("SubspaceIteration :: solveYourselfAt: initial vectors for iteration:\n");
  r.printYourself();
# endif
  
  //ish = 0;
  a-> factorized () ;
  //
  // start of iteration loop
  //
  nite = 0;
  do {                          // label 100
    nite ++;
    ls << "Iteration no. " << nite << endl_notime;
# ifdef DETAILED_REPORT
    printf("SubspaceIteration :: solveYourselfAt: Iteration loop no. %d\n",nite);
# endif
    //
    // compute projection ar and br of matrices a , b
    // 
    {
      ls << "projection...";
      for (j = 1; j<= nc; j++) {
        for (k = 1; k<= nn; k++) tt.at(k) = r.at(k,j) ;
        a->backSubstitutionWith(tt) ;
        for (i = j; i<= nc; i++) {
          art = 0.;
          for (k = 1; k<= nn; k++) art += r.at(k,i)*tt.at(k);
          ar.at(j,i) = art ;
        }
        for (k = 1; k<= nn; k++) r.at(k,j) = tt.at(k);
      }
      ar.symmetrized();            // label 110
#if defined( DETAILED_REPORT ) && DETAILED_REPORT > 1
      printf("SubspaceIteration :: solveYourselfAt: Printing projection matrix ar\n");
      ar.printYourself ();
#endif
      //
      for (j = 1; j<= nc; j++) {   
        for (k = 1; k<= nn ; k++) tt.at(k)=r.at(k,j);
        b->times(tt, temp);
        for (i = j; i<= nc; i++) {
          brt = 0.;
          for (k = 1; k<= nn; k++) brt += r.at(k,i)*temp.at(k) ;
          br.at(j,i) = brt ;
        }                         // label 180
        
        if (iconv < 0) 
          for (k = 1; k<= nn; k++) r.at(k,j)=temp.at(k);
        //delete temp;
      }                           // label 160
      br.symmetrized();
#if defined( DETAILED_REPORT ) && DETAILED_REPORT > 1
      printf("SubspaceIteration :: solveYourselfAt: Printing projection matrix br\n");
      br.printYourself ();
#endif
      ls << "done" << endl;
    }    
    // 
    // solution of reduced eigenvalue problem
    //
    {
      ls << "reduced problem solution...";
      mtd -> solve (&ar,&br,&eigv,&vec);
      ls << "done" << endl;
    }
    //
    // sorting eigenvalues acording to their values
    //
    do {
      is = 0; // label 350
      for (i = 1; i<= nc1; i++) {
        if (fabs(eigv.at(i+1)) < fabs(eigv.at(i))) {
          is++;
          eigvt = eigv.at(i+1);
          eigv.at(i+1) = eigv.at(i) ;
          eigv.at(i)   = eigvt ;
          for (k = 1; k<= nc; k++) {
            rt = vec.at(k,i+1) ;
            vec.at(k,i+1) = vec.at(k,i) ;
            vec.at(k,i)   = rt ;
          }
        }
      }                         // label 360
    } while (is != 0);
# ifdef DETAILED_REPORT
    printf ("SubspaceIteration :: solveYourselfAt: current eigen values of reduced problem \n");
    eigv.printYourself();
# endif
    //
    // compute eigenvectors 
    //
    for (i = 1; i<= nn; i++) { // label 375
      for (j = 1; j<= nc; j++) tt.at(j) = r.at(i,j) ;
      for (k = 1; k<= nc; k++) {
        rt = 0.;
        for (j = 1; j<= nc; j++) rt += tt.at(j)*vec.at(j,k) ;
        r.at(i,k) = rt ;
      }
    }                           // label 420
    //
    // convergency check
    //
    if (iconv >  0) break;
    for (i = 1; i<= nc; i++) {
      dif = (eigv.at(i) - d.at(i));
      rtolv.at(i) = fabs(dif / eigv.at(i));
    }
# ifdef DETAILED_REPORT
    printf ("SubspaceIteration :: solveYourselfAt: Reached precision of eigenvalues:\n");
    rtolv.printYourself();
# endif
    double max_rtolv = 0;
    for (i = 1; i<= nroot; i++) {
      max_rtolv = max(max_rtolv, rtolv.at(i));
    }
    ls << "max eigenvalue error=" << max_rtolv << " (tolerance=" << rtol << ")" << endl_notime;
    if (max_rtolv > rtol) goto label400 ;
    converged = true;
    for (i = 1; i<= nroot; i++) {
    	ls << "		eigenvalue " << i << " error=" << rtolv.at(i) << endl_notime;
    }
    break;
  label400:
    if (nite >= nitem) {
      converged = false;
      for (i = 1; i<= nroot; i++) {
    	ls << "		eigenvalue " << i << " error=" << rtolv.at(i) << endl_notime;
      }
      iconv = 2;
      break;
    }                          
    for (i = 1; i<= nc ; i++) d.at(i) = eigv.at(i); // label 410 and 440
    continue;
  } while (1);
  //
  // May be good idea to compute som norm of erros
  //
  // copy result into caller's arrays
  for (i = 1; i<= nroot; i++) {
    _eigv->at(i) = eigv.at(i);
    for (j=1; j<= nn; j++) _r->at(j,i) = r.at(j,i);
  }
  solved = 1;
  return converged;
}



}
