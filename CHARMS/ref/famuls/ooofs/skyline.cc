#include "skyline.h"
#include "flotmtrx.h"
#include "intarray.h"
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <iostream>
#include "logger_stream_mgr.h"


namespace OOOFS {

  Skyline :: Skyline (int n) : SparseMtrx (n,n) 
  {
    nwk          = 0;
    adr          = NULL;
    mtrx         = NULL;
    isFactorized = false;
  }


  Skyline :: Skyline () : SparseMtrx ()
  {
    nwk          = 0 ;
    adr          = NULL ;
    mtrx         = NULL ;
    isFactorized = false ;
  }


  Skyline :: ~Skyline ()
  {
    // Destructor.
    if (this->giveNumberOfRows()) {
      delete [] mtrx;
      delete (adr);
    }
  }


  double&    
    Skyline ::   at (int i,int j) 
  {
    // returns (i,j) element of the receiver
    // indexes are checked if DEBUG is true

    int d1,k,ind;

#ifdef DEBUG
    // check size
    if ((i > this->giveNumberOfRows()) || (j>this->giveNumberOfRows())){
      printf ("skyline error in function at: dimension mismatch \nfile %s, line %d\n",
        __FILE__,__LINE__) ;
      exit(0) ;}
#endif
    // only upper triangular part of skyline is stored
    if (j<i) {k=i; i=j; j=k;}

    d1 = this->adr->at(i);
    ind = d1 + (j-i);

    if ((adr->at(i+1)-adr->at(i)) <= (j-i)) {
      printf("ERROR: skyline : at - request for element which is not in sparse mtrx\nfile %s, line %d\n",
        __FILE__,__LINE__);
      exit(1);
      //
      // NOTE:
      //
      // don't return reference to some zero value; it is true, but possible change
      // of its value will require rebuilding internal storage structure
      // of sparse matrix
      //
    }

    return mtrx[ind];
  }

  double
    Skyline ::   at (int i,int j) const
  {
    // returns (i,j) element of the receiver
    // indexes are checked if DEBUG is true

    int d1,k,ind;

#ifdef DEBUG
    // check size
    if ((i > this->giveNumberOfRows()) || (j>this->giveNumberOfRows())){
      printf ("skyline error in function at: dimension mismatch \nfile %s, line %d\n",
        __FILE__,__LINE__) ;
      exit(0) ;}
#endif
    // only upper triangular part of skyline is stored
    if (j<i) {k=i; i=j; j=k;}

    d1 = this->adr->at(i);
    ind = d1 + (j-i);

    if ((adr->at(i+1)-adr->at(i)) <= (j-i)) {
      printf("ERROR: skyline : at - request for element which is not in sparse mtrx\nfile %s, line %d\n",
        __FILE__,__LINE__);
      exit(1);
      //
      // NOTE:
      //
      // don't return reference to some zero value; it is true, but possible change
      // of its value will require rebuilding internal storage structure
      // of sparse matrix
      //
    }

    return mtrx[ind];
  }


  void
    Skyline :: toFloatMatrix (FloatMatrix& answer) const
  {
    // Returns a matrix, the receiver in a full storage form. This is useful
    // for debugging and printings.

    //FloatMatrix* answer ;
    int          i,j,d1,d2,pk, size ;
    // double       coeff ;

    size = this->giveNumberOfRows ();
    //answer = new FloatMatrix(size, size) ;
    answer.resize (size, size) ;
    answer.zero();

    for (j=1 ; j<=size ; j++) {
      d1 = adr->at(j);
      d2 = adr->at(j+1);
      pk = j;
      for (i=d1 ; i<d2; i++) {
        answer.at(pk,j)=mtrx[i];
        pk--;
      }
    }
    return  ;
  }



  int  Skyline :: assemble  (const IntArray& loc, const FloatMatrix& mat) 
  {
    // Assembles the elemental matrix 'mat' to the receiver, using 'loc' as a
    // location array. The values in ke corresponding to a zero coefficient
    // in loc are not assembled.

    //  IntArray loc ;
    //  FloatMatrix mat;
    //	int ielem,i,j,ac,ac1,ac2,ndofe;
    //	Domain* domain = eModel->giveDomain();
    int i,j,ac,ac1,ac2,ndofe;

    /*
    int nelem = domain -> giveNumberOfElements ();
    for ( ielem = 1; ielem <= nelem ; ielem++ ) {
    domain -> giveElement(ielem) -> giveLocationArray (loc);
    eModel->giveElementCharacteristicMatrix(mat, ielem, type, tStep );
    */

#  ifdef DEBUG
    int dim = mat.giveNumberOfRows() ;
    if (dim != loc.giveSize()) {
      printf ("error : dimension of 'mat' and 'loc' mismatch \nfile %s, line %d\n",
        __FILE__,__LINE__) ;
      exit(0) ;}
#  endif

    ndofe = mat.giveNumberOfRows() ;

    for (i=1;i<=ndofe;i++){
      ac1=loc.at(i);
      if (ac1==0)  continue;
      for (j=1;j<=ndofe;j++){
        ac2=loc.at(j);
        if (ac2==0) continue;
        if (ac1>ac2) continue;
        ac=adr->at(ac2)+ac2-ac1;
        mtrx[ac]+=mat.at(i,j);
      }
    }

    return 1;
  }




  int  Skyline :: assemble (const IntArray& rloc, const IntArray& cloc, const FloatMatrix& mat) 
  {
    printf ("error : assemble of 'mat' using  'rloc and cloc' unsupported \nfile %s, line %d\n",
      __FILE__,__LINE__) ;
    exit(0) ;
    return 0;
  }


  FloatArray*  Skyline :: backSubstitutionWith (FloatArray& y) const
    // Returns the solution x of the system U.x = y , where U is the receiver.
    // nota : x overwrites y
  {
    // allocation of answer
    FloatArray solution (y.giveSize());
    int i,k,ack,ack1,acs,n;
    int size = this->giveNumberOfRows();
    double s ;

    //solution = new FloatArray (y.giveSize());


    /*****************************/
    /*  modifikace prave strany  */
    /*****************************/
    n = size;
    for (k=2;k<=n;k++){
      ack=adr->at(k);  ack1=adr->at(k+1);
      s=0.0;  acs=k-(ack1-ack)+1;
      for (i=ack1-1;i>ack;i--){
        s+=mtrx[i]*y.at(acs);
        acs++;
      }
      y.at(k)-=s;
    }
    /*****************/
    /*  zpetny chod  */
    /*****************/
    for (k=1;k<=n;k++){
      acs=adr->at(k);
      y.at(k)/=mtrx[acs];
    }

    for (k=n;k>0;k--){
      ack=adr->at(k);  ack1=adr->at(k+1);
      solution.at(k)=y.at(k);
      acs=k-(ack1-ack)+1;
      for (i=ack1-1;i>ack;i--){
        y.at(acs)-=mtrx[i]*solution.at(k);
        acs++;
      }
    }
    y = solution;
    return &y;
  }



  int Skyline :: buildInternalStructure (LOGGER_STREAM &ls, int tot_of_equations, int col_height[]) 
  {
    // first create array of 
    // maximal column height for assembled characteristics matrix
    //

    int j,js,ieq,maxle ;  
    int i,ac1;
    int neq = tot_of_equations;

    IntArray loc;
    IntArray* mht = new IntArray (neq);

    for ( j =1 ; j<=neq; j++) mht->at(j)=col_height[j-1];          // initialize column height

    // increases number of columns according to size of mht
    // mht is array containing minimal equation number per column
    // This method also increases column height.

    if (this->adr) delete adr;
    adr = new IntArray (neq+1);

    ac1 = 1;
    for (i=1;i<=neq;i++){
      adr->at(i) = ac1;
      ac1+= (i-mht->at(i)+1);
    }
    adr->at(neq+1) = ac1;
    nRows = nColumns = neq;
    nwk  = ac1;
    if (mtrx) delete [] mtrx;
    ls << "allocating " << ac1 << " words for skyline matrix" << endl_notime;
    mtrx = new double [ac1];

    delete mht;

    return true;
  }



  SparseMtrx*  Skyline :: factorized ()
  {
    // Returns the receiver in  U(transp).D.U  Crout factorization form.

    int i,j,k,aci,aci1,acj,acj1,ack,ack1,ac,acs,acri,acrk,n;
    double s,g;

    /**********************/
    /*  eliminace matice  */
    /**********************/
    if (isFactorized) return this;
    n = this->giveNumberOfRows();

    for (k=2;k<=n;k++){
      /*  smycka pres sloupce matice  */
      ack=adr->at(k);  ack1=adr->at(k+1);
      acrk=k-(ack1-ack)+1;
      for (i=acrk+1;i<k;i++){
        /*  smycka pres prvky jednoho sloupce matice  */
        aci=adr->at(i);  aci1=adr->at(i+1);
        acri=i-(aci1-aci)+1;
        if (acri<acrk)  ac=acrk;
        else            ac=acri;
        acj=k-ac+ack;  acj1=k-i+ack;
        acs=i-ac+aci;  s=0.0;
        for (j=acj;j>acj1;j--){
          s+=mtrx[j]*mtrx[acs];
          acs--;
        }
        mtrx[acj1]-=s;
      }
      /*  uprava diagonalniho prvku  */
      s=0.0;
      for (i=ack1-1;i>ack;i--){
        g=mtrx[i];
        acs=adr->at(acrk);
        acrk++;
        mtrx[i]/=mtrx[acs];
        s+=mtrx[i]*g;
      }
      mtrx[ack]-=s;
    }

    isFactorized = true ;

    return this ;
  }



  void  Skyline :: times (const FloatArray& x, FloatArray& answer) const
  {
    // Computes y, the results  of the  y = U.x, where U is
    // the receiver. Returns the result.

    int i,j,k,acb,acc,aci,aci1,ac,n;
    double s;

    // 
    // first check sizes
    //
    if (this->giveNumberOfRows() != (n = x.giveSize())) {
      printf ("Class Skyline, function Times : size mismatch\nfile %s, line %d\n",
        __FILE__,__LINE__);
      exit(1);
    }
    //FloatArray *answer = new FloatArray(n);
    answer.resize (n);
    answer.zero();

    acc=1;	for (i=1;i<=n;i++){
      aci=adr->at(i);  aci1=adr->at(i+1);
      ac=i-(aci1-aci)+1;
      s=0.0;  acb=ac;
      for (k=aci1-1;k>=aci;k--){
        s+=mtrx[k]*x.at(acb);
        acb++;
      }
      answer.at(acc)=s;
      acc++;

      for (j=ac;j<i;j++){
        aci1--;
        s=mtrx[aci1];
        answer.at(j)+=s*x.at(i);
        aci++;
      }
    }
    return ;
  }


  void Skyline :: times (double x) 
  {
    // Multiplies receiver by scalar value.

    int j ;

    for (j=0 ; j< nwk ; j++)
      mtrx[j] *= x;

  }


  void  Skyline :: printYourself () const 
  {
    // Prints the receiver on screen.
    FloatMatrix copy ;

    this -> toFloatMatrix(copy) ;
    copy.printYourself() ;
  }


  SparseMtrx*  Skyline :: zero ()
  {
    // Returns the receiver with all coefficients set to zero.

    int j ;

    for (j=0 ; j< nwk ; j++)
      mtrx[j] = 0.0;
    isFactorized = false ;

    return this ;
  }

  SparseMtrx* Skyline :: GiveCopy () const
  {

    Skyline* answer;
    IntArray* adr1;
    double*   mtrx1;
    int neq, i;

    neq = this->giveNumberOfRows();
    adr1 = new IntArray (neq + 1);

    for (i=1; i<=neq+1; i++){
      adr1->at(i) = this->adr->at(i);
    }

    mtrx1 = new double[ (this->nwk)];
    for (i=0; i< this->nwk; i++)
      mtrx1[i]=this->mtrx[i];

    answer = new Skyline (neq, this->nwk, mtrx1, adr1);

    return answer;

  }

  Skyline :: Skyline (int neq, int nwk1, double* mtrx1, IntArray* adr1) : SparseMtrx (neq, neq)
  {
    // constructor
    // sets internal member data to given parameters
    // used only by GiveCopy() member function

    nwk  = nwk1;
    mtrx = mtrx1;
    adr  = adr1;
  }



}



