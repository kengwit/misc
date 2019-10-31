#ifndef BTCBBBAR_H
#  define BTCBBBAR_H

/**
   The expressions here had been generated from

syms N_x N_y N_z real
syms aN_x aN_y aN_z real
syms aN_x_avg aN_y_avg aN_z_avg real
B(1,1)=sym('N_x'); B(1,2)=sym('0'); B(1,3)=sym('0');
B(2,1)=sym('0'); B(2,2)=sym('N_y'); B(2,3)=sym('0'); 
B(3,1)=sym('0'); B(3,2)=sym('0'); B(3,3)=sym('N_z'); 
B(4,1)=sym('N_y'); B(4,2)=sym('N_x'); B(4,3)=sym('0');
B(5,1)=sym('0'); B(5,2)=sym('N_z'); B(5,3)=sym('N_y');
B(6,1)=sym('N_z'); B(6,2)=sym('0'); B(6,3)=sym('N_x');
Bdil=[sym('N_x') sym('N_y') sym('N_z');...
      sym('N_x') sym('N_y') sym('N_z');...
      sym('N_x') sym('N_y') sym('N_z');...
      sym('0') sym('0') sym('0');...
      sym('0') sym('0') sym('0');...
      sym('0') sym('0') sym('0')]/3;
Bbardil=[sym('N_x_avg') sym('N_y_avg') sym('N_z_avg');...
      sym('N_x_avg') sym('N_y_avg') sym('N_z_avg');...
      sym('N_x_avg') sym('N_y_avg') sym('N_z_avg');...
      sym('0') sym('0') sym('0');...
      sym('0') sym('0') sym('0');...
      sym('0') sym('0') sym('0')]/3;
Bdev=B-Bdil;
Bbar=Bdev+Bbardil;

aB(1,1)=sym('aN_x'); aB(1,2)=sym('0'); aB(1,3)=sym('0');
aB(2,1)=sym('0'); aB(2,2)=sym('aN_y'); aB(2,3)=sym('0'); 
aB(3,1)=sym('0'); aB(3,2)=sym('0'); aB(3,3)=sym('aN_z'); 
aB(4,1)=sym('aN_y'); aB(4,2)=sym('aN_x'); aB(4,3)=sym('0');
aB(5,1)=sym('0'); aB(5,2)=sym('aN_z'); aB(5,3)=sym('aN_y');
aB(6,1)=sym('aN_z'); aB(6,2)=sym('0'); aB(6,3)=sym('aN_x');
aBdil=[sym('aN_x') sym('aN_y') sym('aN_z');...
      sym('aN_x') sym('aN_y') sym('aN_z');...
      sym('aN_x') sym('aN_y') sym('aN_z');...
      sym('0') sym('0') sym('0');...
      sym('0') sym('0') sym('0');...
      sym('0') sym('0') sym('0')]/3;
aBbardil=[sym('aN_x_avg') sym('aN_y_avg') sym('aN_z_avg');...
      sym('aN_x_avg') sym('aN_y_avg') sym('aN_z_avg');...
      sym('aN_x_avg') sym('aN_y_avg') sym('aN_z_avg');...
      sym('0') sym('0') sym('0');...
      sym('0') sym('0') sym('0');...
      sym('0') sym('0') sym('0')]/3;
aBdev=aB-aBdil;
aBbar=aBdev+aBbardil;

C_iso(1,1)=sym('lambda + 2*mu'); C_iso(1,2)=sym('lambda');        C_iso(1,3)=sym('lambda');        C_iso(1,4)=sym('0');  C_iso(1,5)=sym('0');  C_iso(1,6)=sym('0');
C_iso(2,1)=sym('lambda');        C_iso(2,2)=sym('lambda + 2*mu'); C_iso(2,3)=sym('lambda');        C_iso(2,4)=sym('0');  C_iso(2,5)=sym('0');  C_iso(2,6)=sym('0');
C_iso(3,1)=sym('lambda');        C_iso(3,2)=sym('lambda');        C_iso(3,3)=sym('lambda + 2*mu'); C_iso(3,4)=sym('0');  C_iso(3,5)=sym('0');  C_iso(3,6)=sym('0');
C_iso(4,1)=sym('0');             C_iso(4,2)=sym('0');             C_iso(4,3)=sym('0');             C_iso(4,4)=sym('mu'); C_iso(4,5)=sym('0');  C_iso(4,6)=sym('0'); 
C_iso(5,1)=sym('0');             C_iso(5,2)=sym('0');             C_iso(5,3)=sym('0');             C_iso(5,4)=sym('0');  C_iso(5,5)=sym('mu'); C_iso(5,6)=sym('0'); 
C_iso(6,1)=sym('0');             C_iso(6,2)=sym('0');             C_iso(6,3)=sym('0');             C_iso(6,4)=sym('0');  C_iso(6,5)=sym('0');  C_iso(6,6)=sym('mu'); 

C(1,1)=sym('C[0][0]');
C(1,2)=sym('C[0][1]');
C(1,3)=sym('C[0][2]');
C(1,4)=sym('C[0][3]');
C(1,5)=sym('C[0][4]');
C(1,6)=sym('C[0][5]');
C(2,1)=sym('C[1][0]');
C(2,2)=sym('C[1][1]');
C(2,3)=sym('C[1][2]');
C(2,4)=sym('C[1][3]');
C(2,5)=sym('C[1][4]');
C(2,6)=sym('C[1][5]');
C(3,1)=sym('C[2][0]');
C(3,2)=sym('C[2][1]');
C(3,3)=sym('C[2][2]');
C(3,4)=sym('C[2][3]');
C(3,5)=sym('C[2][4]');
C(3,6)=sym('C[2][5]');
C(4,1)=sym('C[3][0]');
C(4,2)=sym('C[3][1]');
C(4,3)=sym('C[3][2]');
C(4,4)=sym('C[3][3]');
C(4,5)=sym('C[3][4]');
C(4,6)=sym('C[3][5]');
C(5,1)=sym('C[4][0]');
C(5,2)=sym('C[4][1]');
C(5,3)=sym('C[4][2]');
C(5,4)=sym('C[4][3]');
C(5,5)=sym('C[4][4]');
C(5,6)=sym('C[4][5]');
C(6,1)=sym('C[5][0]');
C(6,2)=sym('C[5][1]');
C(6,3)=sym('C[5][2]');
C(6,4)=sym('C[5][3]');
C(6,5)=sym('C[5][4]');
C(6,6)=sym('C[5][5]');
*/

static void
btcb_3d_aniso_symm_C_bbar (double aN_x, double aN_y, double aN_z,
			   double N_x, double N_y, double N_z,                      
			   double N_x_avg, double N_y_avg, double N_z_avg,                      
			   double C[6][6], double mult, double K[3][3])
{
  double CBbar[6][3];
  double Bbar[6][3];
  Bbar[0][0] = 2.0/3.0*N_x+N_x_avg/3.0;      Bbar[0][1] = -N_y/3.0+N_y_avg/3.0; Bbar[0][2] = -N_z/3.0+N_z_avg/3.0;
  Bbar[1][0] = -N_x/3.0+N_x_avg/3.0;    Bbar[1][1] = 2.0/3.0*N_y+N_y_avg/3.0;  Bbar[1][2] = -N_z/3.0+N_z_avg/3.0; 
  Bbar[2][0] = -N_x/3.0+N_x_avg/3.0;  Bbar[2][1] = -N_y/3.0+N_y_avg/3.0;Bbar[2][2] = 2.0/3.0*N_z+N_z_avg/3.0; 
  Bbar[3][0] = N_y;    Bbar[3][1] = N_x;   Bbar[3][2] = 0.0; 
  Bbar[4][0] = 0.0;      Bbar[4][1] = N_z;  Bbar[4][2] = N_y; 
  Bbar[5][0] = N_z;  Bbar[5][1] = 0.0;  Bbar[5][2] = N_x;
  double aBbar[6][3];
  aBbar[0][0] = 2.0/3.0*aN_x+N_x_avg/3.0;     aBbar[0][1] = -aN_y/3.0+N_y_avg/3.0;   aBbar[0][2] = -aN_z/3.0+N_z_avg/3.0;  
  aBbar[1][0] = -aN_x/3.0+N_x_avg/3.0; aBbar[1][1] = 2.0/3.0*aN_y+N_y_avg/3.0; aBbar[1][2] = -aN_z/3.0+N_z_avg/3.0; 
  aBbar[2][0] = -aN_x/3.0+N_x_avg/3.0;  aBbar[2][1] = -aN_y/3.0+N_y_avg/3.0;      aBbar[2][2] = 2.0/3.0*aN_z+N_z_avg/3.0;      
  aBbar[3][0] = aN_y;      aBbar[3][1] = aN_x;      aBbar[3][2] = 0.0;    
  aBbar[4][0] = 0.0;      aBbar[4][1] = aN_z;      aBbar[4][2] = aN_y; 
  aBbar[5][0] = aN_z;      aBbar[5][1] = 0.0;      aBbar[5][2] = aN_x;
  for (int i=0; i<6; i++) {
    for (int j=0; j<3; j++) {
      double sum=0;
      for (int k=0; k<6; k++) {
	sum += C[i][k]*Bbar[k][j];
      }
      CBbar[i][j]=sum;
    }
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      double sum=0;
      for (int k=0; k<6; k++) {
	sum += aBbar[k][i]*CBbar[k][j];
      }
      K[i][j]=mult * sum;
    }
  }
}

static void
btcb_3d_iso_symm_C_bbar (double aN_x, double aN_y, double aN_z,
			 double N_x, double N_y, double N_z,
			 double N_x_avg, double N_y_avg, double N_z_avg,                      
			 double lambda, double mu, double mult, double K[3][3])
{
  /* Produced by ccode(simplify(aBbar'*C_iso*Bbar)) */
  K[0][0] = mult * (4.0/3.0*aN_x*mu*N_x+N_x_avg*lambda*N_x_avg+2.0/3.0*N_x_avg*mu*N_x_avg+aN_y*mu*N_y+aN_z*mu*N_z);      
  K[0][1] = mult * (-2.0/3.0*aN_x*mu*N_y+N_x_avg*lambda*N_y_avg+2.0/3.0*N_x_avg*mu*N_y_avg+aN_y*mu*N_x);      
  K[0][2] = mult * (-2.0/3.0*aN_x*mu*N_z+N_x_avg*lambda*N_z_avg+2.0/3.0*N_x_avg*mu*N_z_avg+aN_z*mu*N_x);
  K[1][0] = mult * (-2.0/3.0*aN_y*mu*N_x+N_y_avg*lambda*N_x_avg+2.0/3.0*N_y_avg*mu*N_x_avg+aN_x*mu*N_y); 
  K[1][1] = mult * (4.0/3.0*aN_y*mu*N_y+N_y_avg*lambda*N_y_avg+2.0/3.0*N_y_avg*mu*N_y_avg+aN_x*mu*N_x+aN_z*mu*N_z); 
  K[1][2] = mult * (-2.0/3.0*aN_y*mu*N_z+N_y_avg*lambda*N_z_avg+2.0/3.0*N_y_avg*mu*N_z_avg+aN_z*mu*N_y); 
  K[2][0] = mult * (-2.0/3.0*aN_z*mu*N_x+N_z_avg*lambda*N_x_avg+2.0/3.0*N_z_avg*mu*N_x_avg+aN_x*mu*N_z);  
  K[2][1] = mult * (-2.0/3.0*aN_z*mu*N_y+N_z_avg*lambda*N_y_avg+2.0/3.0*N_z_avg*mu*N_y_avg+aN_y*mu*N_z);   
  K[2][2] = mult * (4.0/3.0*aN_z*mu*N_z+N_z_avg*lambda*N_z_avg+2.0/3.0*N_z_avg*mu*N_z_avg+aN_y*mu*N_y+aN_x*mu*N_x);
}

/**
   Transpose K in place.
*/
static void
transp (double K[3][3])
{
#define SWP(a,b) { double c = a; a = b; b = c; }
  SWP(K[0][1], K[1][0]);
  SWP(K[0][2], K[2][0]);
  SWP(K[1][2], K[2][1]);
}

#endif
