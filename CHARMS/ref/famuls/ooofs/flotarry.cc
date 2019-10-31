#include "flotarry.h"
#include "intarray.h"
#include "flotmtrx.h"
#include <math.h>
#include <string.h>

#if !defined(max)
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#if !defined(min)
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

namespace OOOFS {

  FloatArray :: FloatArray (int n)
   // Constructor : creates an array of size n (filled with garbage).
{
   allocatedSize = size = n ;
   if (size)
      values = new double [size];
   else
      values = NULL ;
}

FloatArray :: FloatArray (const FloatArray & src) 
{
// copy constructor

	double * srcVal;

	allocatedSize = size = src.size ;
	if (size)
		values = new double[size] ;
	else
		values = NULL ;
	
	srcVal = src.givePointer();
	for (int i=0; i< size; i++) this->values[i]=srcVal[i];
}

FloatArray&
FloatArray :: operator= (const FloatArray& src)
{
	// assignment: cleanup and copy
	double *srcVal;

	if (this != &src) { // beware of s=s;
		this->resize (src.size);
		
		srcVal = src.givePointer();
		for (int i=0; i< size; i++) this->values[i]=srcVal[i];
	}
	return *this;
}


void
FloatArray :: add (const FloatArray& b)
   // Performs the operation a=a+b, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of b. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (b.giveSize()==0)
      return  ;

   if (! size)                                // null-sized array
		 {resize (b.size); zero();}

	 if (size != b.size) {                  // unmatching sizes
		 printf ("FloatArray dimension mismatch in a[%d]->add(b[%d])\n",
						 size,b.size) ;
		 exit(0) ;}

   p1 = values ;
   p2 = b.values ;
   i  = size ;
   while (i--)
		 *p1++ += *p2++ ;
   return  ;
}

void  FloatArray :: substract (const FloatArray& src)
   // Performs the operation a=a-src, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of src. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (src.giveSize()==0)
      return  ;

   if (! size) {                              // null-sized array
		 this->resize (size = src.size) ;
		 this->zero();
	 }

#  ifdef DEBUG
	 if (size != src.size) {                  // unmatching sizes
		 printf ("FloatArray dimension mismatch in a[%d]->add(b[%d])\n",
						 size,src.size) ;
		 exit(0) ;}
#  endif

   p1 = values ;
   p2 = src.values ;
   i  = size ;
   while (i--)
      *p1++ -= *p2++ ;
   return  ;
}



double  dotProduct (const FloatArray& p1, const FloatArray& p2, register int i)
   // A non-member function. Returns the dot product of the first 'i' coef-
   // ficienst of the two arrays P1 and P2. This method applies to many
   // situations, eg row*column products with matrices.
{
   double answer ;

   answer = 0. ;
   while (i) {
		 answer += p1.at(i) * p2.at(i) ;
		 i--;
	 }
   return answer ;
}

void
FloatArray :: beSubArrayOf (const FloatArray& src, const IntArray &indx)
//
// returns subVector from receiver
// subVector size will be indx max val size
// and on i-th position of subVector will be this->at(indx->at(i))
//
{
	//FloatArray *answer;
	int i,ii,n,isize;

	n = indx.giveSize();

	if (src.size != n) {printf ("beSubArrayOf - size mismatch");
						 exit(0);}

	for (isize=0,i=1; i<=n;i++) if(indx.at(i)>isize) isize = indx.at(i);
	//answer = new FloatArray (isize);
	this->resize (isize);
	for (i=1; i<= n; i++) {
		ii = indx.at(i);
		if (ii > 0) this->at(ii) = src.at(i);
	}
	return ;
}


void 
FloatArray :: addSubVector (const FloatArray& src, int si)
{
	int i, reqSize, n = src.giveSize();

	si --;
	reqSize = si+n;
	if (this->giveSize() < reqSize) this->resize(reqSize);

	for (i=1; i<= n; i++) this->at(si+i) += src.at(i);
}



void
FloatArray :: beVectorProductOf (const FloatArray& v1, const FloatArray& v2)
//
// computes vector product v1 x v2
// and stores result into receiver
//
{
	// check proper bunds
	if ((v1.giveSize() != 3) || (v2.giveSize() != 3)) 
		{
			printf ("error in FloatArray::VectorProduct size mismatch, size is not equal to 3 \n") ;
			exit(0) ;
		}
	this->resize (3);
	
	this->at(1) = v1.at(2)*v2.at(3)-v1.at(3)*v2.at(2);
	this->at(2) = v1.at(3)*v2.at(1)-v1.at(1)*v2.at(3);
	this->at(3) = v1.at(1)*v2.at(2)-v1.at(2)*v2.at(1);
	
	return ;
}

double
FloatArray :: distance (const FloatArray &from) const
//
// returns distance between receiver and from from
// computed using generalized pythagora formulae
//
{
	double dist = 0.;

	if (size != from.giveSize())  {
		printf ("error in FloatArray::distance size mismatch (%d) is not equal to (%d) \n",
				  size,from.giveSize()) ;
		exit(0) ;}
	
	for (int i = 1; i<= size; i++)
		dist += (this->at(i)-from.at(i))*(this->at(i)-from.at(i));
	
	return sqrt(dist);
	
} 


/*******************************************************************************/

FloatArray*  FloatArray :: add (FloatArray* b)
   // Performs the operation a=a+b, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of b. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (!b || b->giveSize()==0)
      return this ;

   if (! size) {                              // null-sized array
      size   = b -> size ;
      values = new double[(size)] ;}

#  ifdef DEBUG
      if (size != b->size) {                  // unmatching sizes
	 printf ("FloatArray dimension mismatch in a[%d]->add(b[%d])\n",
		  size,b->size) ;
	 exit(0) ;}
#  endif

   p1 = values ;
   p2 = b -> values ;
   i  = size ;
   while (i--)
      *p1++ += *p2++ ;
   return this ;
}


FloatArray*  FloatArray :: substract (FloatArray* b)
   // Performs the operation a=a-b, where a stands for the receiver. If the
   // receiver's size is 0, adjusts its size to that of b. Returns the
   // receiver.
{
   register int i ;
   double       *p1,*p2 ;

   if (!b || b->giveSize()==0)
      return this ;

   if (! size) {                              // null-sized array
      size   = b -> size ;
      values = new double[(size)] ;}

#  ifdef DEBUG
      if (size != b->size) {                  // unmatching sizes
	 printf ("FloatArray dimension mismatch in a[%d]->add(b[%d])\n",
		  size,b->size) ;
	 exit(0) ;}
#  endif

   p1 = values ;
   p2 = b -> values ;
   i  = size ;
   while (i--)
      *p1++ -= *p2++ ;
   return this ;
}



void  FloatArray :: assemble (const FloatArray& fe, const IntArray& loc)
   // Assembles the array fe (typically, the load vector of a finite
   // element) to the receiver, using loc as location array.
{
   int i,ii,n ;

#  ifdef DEBUG
      if ((n=fe.giveSize()) != loc.giveSize()) {
	 printf ("dimensions of 'fe' and 'loc' mismatch \n") ;
	 exit(0) ; }
      this -> checkSizeTowards(loc) ;
#  endif

   n = fe.giveSize() ;
   for (i=1 ; i<=n ; i++) {
      ii = loc.at(i) ;
      if (ii)                                  // if non 0 coefficient,
	 this->at(ii) += fe.at(i) ; }         // then assemble
}


#ifdef DEBUG
double&  FloatArray :: at (int i)
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}

double  FloatArray :: at (int i) const
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}
#endif


#ifdef DEBUG
void  FloatArray :: checkBounds(int i) const
   // Checks that the receiver's size is not smaller than 'i'.
{
   if (i<=0) {
      printf ("array error on index : %d <= 0 \n",i) ;
      exit(0) ; }

   if (i>size) {
      printf ("array error on index : %d > %d \n",i,size) ;
      exit(0) ; }
}
#endif

FloatArray*
FloatArray :: GiveSubArray (IntArray *indx)
//
// returns subVector from receiver
// subVector size will be indx max val size
// and on i-th position of subVector will be this->at(indx->at(i))
//
{
	FloatArray *answer;
	int i,ii,n,isize;

	n = indx->giveSize();

	if (size != n) {printf ("GiveSubArray - size mismatch");
						 exit(0);}

	for (isize=0,i=1; i<=n;i++) if(indx->at(i)>isize) isize = indx->at(i);
	answer = new FloatArray (isize);
	for (i=1; i<= n; i++) {
		ii = indx->at(i);
		if (ii > 0)answer->at(ii) = this->at(i);
	}
	return answer;
}

void  FloatArray :: checkSizeTowards (const IntArray& loc)
   // Expands the receiver if loc points to coefficients beyond the size of
   // the receiver.
{
   int i,n,high ;

   high = 0 ;
   n    = loc.giveSize() ;
   for (i=1 ; i<=n ; i++)
      high = max (high,(loc.at(i))) ;
   if (high > size)                             // receiver must be expanded
      this -> resize (high) ;
}


void  FloatArray :: resize (int n, int allocChunk)
   // Expands the receiver up to size n (n is assumed larger than 'size').
   // Initializes all new coefficients to zero.
{
   register int i ;
   double       *newValues,*p1,*p2 ;

   if (n <= allocatedSize) {size = n; return;}

	 if (allocChunk < 0) allocChunk = 0;
   newValues = new double[(n + allocChunk)] ;

   p1 = values ;
   p2 = newValues ;
   i  = size ;
   while (i--)
      *p2++ = *p1++ ;

   if (values)
		 delete [] (values) ;
   values = newValues ;
   allocatedSize = n + allocChunk;
	 size   = n ;
 }


void  FloatArray :: hardResize (int n)
   // Realocates the receiver with new size.
   // Initializes all new coefficients to zero.
{
   register int i ;
   double       *newValues,*p1,*p2 ;

   // if (n <= allocatedSize) {size = n; return;}

   newValues = new double[(n)] ;

   p1 = values ;
   p2 = newValues ;
   i  = min(size, n) ;
   while (i--)
      *p2++ = *p1++ ;

   if (values)
      delete [] (values) ;
   values = newValues ;
   allocatedSize = size   = n ;
}



int  FloatArray :: containsOnlyZeroes () const
   // Returns True if all coefficients of the receiver are 0, else returns
   // False.
{
   register int i ;
   double       *p ;

   p = values ;
   i = size ;
   while (i--)
      if (*p++ != 0.)
        return false ;

   return true ;
}



void FloatArray :: zero ()
// zeroes all values to zero
{
  int i;
  if (values)
	 for (i=0; i< size; i++) values[i]=0.;
  // return *this ;
}


void  FloatArray :: beProductOf (const FloatMatrix& aMatrix, const FloatArray& anArray)
   // Stores the product of aMatrix * anArray in to receiver
{
   int         i,j, nColumns, nRows ;
   double      sum ;
   //FloatArray* answer ;

#  ifdef DEBUG
      if ((nColumns=aMatrix.giveNumberOfColumns()) - anArray.giveSize()) {
	 printf ("error in product A*x : dimension mismatch \nfile %s, line %d\n",
					 __FILE__,__LINE__) ;
	 this    -> printYourself() ;
	 anArray.printYourself() ;
	 exit(0) ;}
#  endif

	 nColumns = aMatrix.giveNumberOfColumns();
   this->resize (nRows = aMatrix.giveNumberOfRows()) ;
   for (i=1 ; i<=nRows ; i++) {
      sum = 0. ;
      for (j=1 ; j<=nColumns ; j++)
				sum += aMatrix.at(i,j) * anArray.at(j) ;
      this->at(i) = sum ;}
   return ;
}

void  FloatArray :: beTProductOf (const FloatMatrix& aMatrix, const FloatArray& anArray)
   // Stores the product of aMatrix^T * anArray in to receiver
{
   int         i,j, nColumns, nRows ;
   double      sum ;
   //FloatArray* answer ;

#  ifdef DEBUG
      if ((nColumns=aMatrix.giveNumberOfRows()) - anArray.giveSize()) {
	 printf ("error in product A*x : dimension mismatch \nfile %s, line %d\n",
					 __FILE__,__LINE__) ;
	 this    -> printYourself() ;
	 anArray.printYourself() ;
	 exit(0) ;}
#  endif

	 nColumns = aMatrix.giveNumberOfRows();
   this->resize (nRows = aMatrix.giveNumberOfColumns()) ;
   for (i=1 ; i<=nRows ; i++) {
      sum = 0. ;
      for (j=1 ; j<=nColumns ; j++)
				sum += aMatrix.at(j,i) * anArray.at(j) ;
      this->at(i) = sum ;}
   return ;
}


FloatArray*  FloatArray :: negated ()
   // Switches the sign of every coefficient of the receiver. Returns the
   // receiver.
{
   register int i ;
   double       x ;
   double*      p ;

   i = size ;
   p = values ;
   while (i--) {
      x    = - *p ;
      *p++ = x ;}
   return this ;
}


void  FloatArray :: printYourself () const
   // Prints the receiver on screen.
{
   printf ("FloatArray of size : %d \n",size) ;
   for (int i=1 ; i<=size ; ++i)
      printf ("%10.3e  ",this->at(i)) ;
   printf ("\n") ;
}

FloatArray*  FloatArray :: setValuesToZero ()
// zeroes all values to zero
{
  int i;
  if (values)
	 for (i=0; i< size; i++) values[i]=0.;
  return this ;
}


FloatArray*  FloatArray :: beCopyOf (FloatArray*arry)
// Returns the receiver initialized according to array arry
// if arry is NULL nothing done
{
  int i;
  double *toVal;
  if (arry) {
	 if (size != arry->giveSize()) {
		delete [] (values);
  		size = arry->giveSize();
		values = new double[(size)] ;
	 }
	 
	 toVal = arry->givePointer();
	 for (i=0; i< size; i++) values[i]=toVal[i];
	 return this ;
  } else return this;
}
  


void
FloatArray :: rotatedWith (FloatMatrix& r, char mode)
   // Returns the receiver 'a' rotated according the change-of-base matrix r.
   // If mode = 't', the method performs the operation  a = r(transp) * a .
   // If mode = 'n', the method performs the operation  a = r * a .
{
   FloatMatrix  rot ;
   FloatArray   rta ;

   if (mode == 't') {
      rot.beTranspositionOf(r) ;
			rta.beProductOf (rot, *this);
		} else if (mode == 'n') {
			rta.beProductOf (r, *this);
		} else {
			fprintf (stderr, "FloatArray :: rotatedWith: unsupported mode");
			exit (1);
		}

	 *this = rta;
}


FloatArray*  FloatArray :: rotatedWith (FloatMatrix* r, char mode)
   // Returns the receiver 'a' rotated according the change-of-base matrix r.
   // If mode = 't', the method performs the operation  a = r(transp) * a .
   // If mode = 'n', the method performs the operation  a = r * a .
{
   // double       *p1,*p2 ;
   FloatMatrix  *rot ;
   FloatArray   *rta ;

   if (mode == 't')
		 rot = r -> GiveTransposition() ;
   else if (mode == 'n') 
		 rot = r ;
	 else {
			fprintf (stderr, "FloatArray :: rotatedWith: unsupported mode");
			exit (1);
	 }

   rta = rot -> Times(this) ;

/*
   p1 = values ;
   p2 = rta -> values ;

   i  = size ;
   while (i--)
      *p1++ = *p2++ ;
*/

	 *this = *rta;

   if (mode == 't')
      delete rot ;
   delete rta ;
   return this ;
}


FloatArray*  FloatArray :: times (double factor)
   // Multiplies every coefficient of the receiver by factor. Answers the
   // modified receiver.
{
   register int i ;
   double*      p ;

   p = values ;
   i = size ;
   while (i--)
      *(p++) *= factor ;
   return this ;
}


FloatArray*  FloatArray :: Times (double factor) const
   // Returns a new array, whose components are those of the receicer, times
   // factor.
{
   register int i ;
   double       *p1,*p2 ;
   FloatArray*  answer ;

   answer = new FloatArray(size) ;
   p1     = values ;
   p2     = answer -> values ;
   i      = size ;
   while (i--)
      *p2++ = factor * (*p1++) ;
   return answer ;
}


double  dotProduct (double* P1, double* P2, register int i)
   // A non-member function. Returns the dot product of the first 'i' coef-
   // ficienst of the two arrays P1 and P2. This method applies to many
   // situations, eg row*column products with matrices.
{
   double answer ;

   answer = 0. ;
   while (i--)
      answer += *P1++ * *P2++ ;
   return answer ;
}


double
FloatArray :: distance (FloatArray *from)
//
// returns distance between receiver and from from
// computed using generalized pythagora formulae
//
{
	double dist = 0.;

	if (size != from-> giveSize())  {
		printf ("error in FloatArray::distance size mismatch (%d) is not equal to (%d) \n",
				  size,from->giveSize()) ;
		exit(0) ;}
	
	for (int i = 1; i<= size; i++)
		dist += (this->at(i)-from->at(i)) * (this->at(i)-from->at(i));
	
	return sqrt(dist);
	
} 

FloatArray*
FloatArray :: VectorProduct (FloatArray* v2)
//
// computes vector product this x v2
// and return newly allocated result
//
{
	// check proper bunds
	if ((size != 3) || (v2->giveSize() != 3)) 
		{
			printf ("error in FloatArray::VectorProduct size mismatch, size is not equal to 3 \n") ;
			exit(0) ;
		}
	FloatArray *answer = new FloatArray(3);
	
	answer->at(1) = this->at(2)*v2->at(3)-this->at(3)*v2->at(2);
	answer->at(2) = this->at(3)*v2->at(1)-this->at(1)*v2->at(3);
	answer->at(3) = this->at(1)*v2->at(2)-this->at(2)*v2->at(1);

	return answer;
}


FloatArray*
FloatArray :: normalize ()
//
// normalizes receiver to have norm equal to 1.0
//
{
	int i;
	double norm = 0.;
	for (i=1; i<= size; i++)
		norm += this->at(i)*this->at(i);
	
	norm = sqrt (norm);
	if (norm < 1.e-80) {
		printf ("error in FloatArray::normalize cannot norm receiver, norm is too small");
		exit(0) ;
	}
	
	this->times(1./norm);
	
	return this;
}





}
