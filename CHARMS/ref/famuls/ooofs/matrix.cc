#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

namespace OOOFS {

void  Matrix :: checkBounds (int i,int j) const
   // Checks that the receiver includes a position (i,j).
{
   if (i<=0) {
      printf ("matrix error on rows : %d < 0 \n",i) ;
      exit(0) ;}
   if (j<=0) {
      printf ("matrix error on columns : %d < 0 \n",j) ;
      exit(0) ;}
   if (i>nRows){
      printf ("matrix error on rows : %d > %d \n",i,nRows) ;
      exit(0) ;}
   if (j>nColumns){
      printf ("matrix error on columns : %d > %d \n",j,nColumns) ;
      exit(0) ;}
}


}
