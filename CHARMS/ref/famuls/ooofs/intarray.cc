#include "intarray.h"
#include <stdio.h>


namespace OOOFS {


IntArray :: IntArray (int n)
   // Constructor : creates an array of size n (filled with garbage).
{
   allocatedSize = size = n ;
   if (n)
      values = new int[(size)] ;
   else
      values = NULL ;
}


IntArray::IntArray (const IntArray& src)
{
	// copy constructor
	allocatedSize = size = src.size ;
	if (size) values = new int [(size)]; else values = NULL;
	
	int *p2        = src.values ;
  int *p3        = values ;

	int i = size ;
	while (i--)
		*p3++ = *p2++ ;

}


IntArray& 
IntArray :: operator=  (const IntArray& src)
{
	// assignment: cleanup and copy
	if (values) delete [] (values);

	
	allocatedSize = size = src.size ;
	if (size) {
		values = new int [(size)] ;
	
		int *p2        = src.values ;
		int *p3        = values ;
	
		int i = size ;
		while (i--)
			*p3++ = *p2++ ;
	} else { values = NULL;}
	return *this;
}


void
IntArray :: zero () 
{
  int *p1        = values ;
	
	int i = size ;
	while (i--)
		*p1++ = 0;
}


#ifdef DEBUG
int&  IntArray :: at (int i)
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}

int  IntArray :: at (int i) const
   // Returns the i-th coefficient of the receiver. Slow but safe.
{
   this -> checkBounds(i) ;
   return values[i-1] ;
}
#endif


#ifdef DEBUG
void  IntArray :: checkBounds (int i) const
   // Checks that the receiver includes an index i.
{
   if (i<0) {
      printf ("array error on index : %d < 0 \n",i) ;
      exit(0) ;}
   if (i>size) {
      printf ("array error on index : %d > %d \n",i,size) ;
      exit(0) ;}
}
#endif

void 
IntArray :: resize (int n, int allocChunk) 
{
	int *p1, *p2, *newValues, i;

	if (n <= allocatedSize) {size = n; return;}

	if (allocChunk < 0) allocChunk = 0;
	newValues = new int [(n + allocChunk)] ;

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


void  IntArray :: followedBy (const IntArray& b, int allocChunk)
		 // Appends the array 'b' the receiver. Returns the receiver.
{
	register int i ;
	int          newSize ;
	int          *newValues,*p1,*p2,*p3 ;
	
	newSize = size + b.size ;
	if (newSize == size)
		return ;

	if (allocChunk < 0) allocChunk = 0;

	if (newSize > allocatedSize) {
		newValues = new int [(newSize+allocChunk)] ;
		p3        = newValues ;

		p1        = values ;
		p2        = b.values ;
		p3        = newValues ;
	
		i = size ;
		while (i--)
			*p3++ = *p1++ ;
		
		i = b.size ;
		while (i--)
			*p3++ = *p2++ ;
		
		if (values)
                  delete [] (values) ;
		values = newValues ;
		allocatedSize = newSize+allocChunk;
		size   = newSize ;
		
	} else {
		
		p1        = values+size ;
		p2        = b.values ;
		
		i = b.size ;
		while (i--)
			*p1++ = *p2++ ;
		
		size   = newSize;
	}
}



void  IntArray :: followedBy (const int b, int allocChunk)
		 // Appends the array 'b' the receiver. Returns the receiver.
{
	register int i ;
	int          newSize ;
	int          *newValues,*p1,*p3 ;
	
	newSize = size + 1 ;

	if (newSize > allocatedSize) {
		newValues = new int [(newSize+allocChunk)] ;

		p1        = values ;
		p3        = newValues ;
	
		i = size ;
		while (i--)
			*p3++ = *p1++ ;
		
		  *p3 = b ;
		
		if (values)
			delete [] (values) ;
		values = newValues ;
		allocatedSize = newSize + allocChunk;
		size   = newSize ;
		
	} else {
		
		*(values+size) = b ;
		size   = newSize;
	}
}

void  IntArray :: printYourself () const
   // Prints the receiver on screen.
{
   printf ("IntArray of size : %d\n",size) ;
   for (int i=1 ; i<=size ; ++i) {
      if (i>15) {
	 printf ("   (other components not printed)") ;
	 break ;}
      else
	 printf ("%d  ",this->at(i)) ;}
   printf ("\n") ;
}




int
IntArray :: findFirstIndexOf (int value)   const
{
	// finds index of value in receiver
	// if such value  does not exists, returns zero index
	int i;
  for (i=0 ; i<size ; i++) {
		if (values[i] == value) return i+1;
	}
	// nothing found
	return 0;
}


}
