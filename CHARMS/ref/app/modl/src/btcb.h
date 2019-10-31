#ifndef BTCB_H
#  define BTCB_H

/**
   The expressions here had been generated from

B(1,1)=sym('N_x'); B(1,2)=sym('0'); B(1,3)=sym('0');
B(2,1)=sym('0'); B(2,2)=sym('N_y'); B(2,3)=sym('0'); 
B(3,1)=sym('0'); B(3,2)=sym('0'); B(3,3)=sym('N_z'); 
B(4,1)=sym('N_y'); B(4,2)=sym('N_x'); B(4,3)=sym('0');
B(5,1)=sym('0'); B(5,2)=sym('N_z'); B(5,3)=sym('N_y');
B(6,1)=sym('N_z'); B(6,2)=sym('0'); B(6,3)=sym('N_x');

aB(1,1)=sym('aN_x'); aB(1,2)=sym('0'); aB(1,3)=sym('0');
aB(2,1)=sym('0'); aB(2,2)=sym('aN_y'); aB(2,3)=sym('0'); 
aB(3,1)=sym('0'); aB(3,2)=sym('0'); aB(3,3)=sym('aN_z'); 
aB(4,1)=sym('aN_y'); aB(4,2)=sym('aN_x'); aB(4,3)=sym('0');
aB(5,1)=sym('0'); aB(5,2)=sym('aN_z'); aB(5,3)=sym('aN_y');
aB(6,1)=sym('aN_z'); aB(6,2)=sym('0'); aB(6,3)=sym('aN_x');

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
btcb_3d_aniso_symm_C (double aN_x, double aN_y, double aN_z,
                      double N_x, double N_y, double N_z,                      
                      double C[6][6], double mult, double K[3][3])
{
  double C00, C01, C02, C03, C04, C05;
  double C10, C11, C12, C13, C14, C15;
  double C20, C21, C22, C23, C24, C25;
  double C30, C31, C32, C33, C34, C35;
  double C40, C41, C42, C43, C44, C45;
  double C50, C51, C52, C53, C54, C55;
  C00=C[0][0];
  C01=C10=C[0][1];
  C02=C20=C[0][2];
  C03=C30=C[0][3];
  C04=C40=C[0][4];
  C05=C50=C[0][5];
  C11=C[1][1];
  C12=C21=C[1][2];
  C13=C31=C[1][3];
  C14=C41=C[1][4];
  C15=C51=C[1][5];
  C22=C[2][2];
  C23=C32=C[2][3];
  C24=C42=C[2][4];
  C25=C52=C[2][5];
  C33=C[3][3];
  C34=C43=C[3][4];
  C35=C53=C[3][5];
  C44=C[4][4];
  C45=C54=C[4][5];
  C55=C[5][5];
  /* Produced by ccode(simplify(aB'*C*B)) */
  K[0][0] = mult * (N_x*aN_x*C00+N_x*aN_y*C30+N_x*aN_z*C50+N_y*aN_x*C03+N_y*aN_y*C33+N_y*aN_z*C53+N_z*aN_x*C05+N_z*aN_y*C35+N_z*aN_z*C55);
  K[0][1] = mult * (N_y*aN_x*C01+N_y*aN_y*C31+N_y*aN_z*C51+N_x*aN_x*C03+N_x*aN_y*C33+N_x*aN_z*C53+N_z*aN_x*C04+N_z*aN_y*C34+N_z*aN_z*C54);
  K[0][2] = mult * (N_z*aN_x*C02+N_z*aN_y*C32+N_z*aN_z*C52+N_y*aN_x*C04+N_y*aN_y*C34+N_y*aN_z*C54+N_x*aN_x*C05+N_x*aN_y*C35+N_x*aN_z*C55);
  K[1][0] = mult * (N_x*aN_y*C10+N_x*aN_x*C30+N_x*aN_z*C40+N_y*aN_y*C13+N_y*aN_x*C33+N_y*aN_z*C43+N_z*aN_y*C15+N_z*aN_x*C35+N_z*aN_z*C45);
  K[1][1] = mult * (N_y*aN_y*C11+N_y*aN_x*C31+N_y*aN_z*C41+N_x*aN_y*C13+N_x*aN_x*C33+N_x*aN_z*C43+N_z*aN_y*C14+N_z*aN_x*C34+N_z*aN_z*C44);
  K[1][2] = mult * (N_z*aN_y*C12+N_z*aN_x*C32+N_z*aN_z*C42+N_y*aN_y*C14+N_y*aN_x*C34+N_y*aN_z*C44+N_x*aN_y*C15+N_x*aN_x*C35+N_x*aN_z*C45);
  K[2][0] = mult * (N_x*aN_z*C20+N_x*aN_y*C40+N_x*aN_x*C50+N_y*aN_z*C23+N_y*aN_y*C43+N_y*aN_x*C53+N_z*aN_z*C25+N_z*aN_y*C45+N_z*aN_x*C55);
  K[2][1] = mult * (N_y*aN_z*C21+N_y*aN_y*C41+N_y*aN_x*C51+N_x*aN_z*C23+N_x*aN_y*C43+N_x*aN_x*C53+N_z*aN_z*C24+N_z*aN_y*C44+N_z*aN_x*C54);
  K[2][2] = mult * (N_z*aN_z*C22+N_z*aN_y*C42+N_z*aN_x*C52+N_y*aN_z*C24+N_y*aN_y*C44+N_y*aN_x*C54+N_x*aN_z*C25+N_x*aN_y*C45+N_x*aN_x*C55);
}

static void
btcb_3d_iso_symm_C (double aN_x, double aN_y, double aN_z,
                    double N_x, double N_y, double N_z,
                    double lambda, double mu, double mult, double K[3][3])
{
  /* Produced by ccode(simplify(aB'*C_iso*B)) */
  K[0][0] = mult * (aN_x*N_x*lambda+2.0*aN_x*mu*N_x+aN_y*mu*N_y+aN_z*mu*N_z);
  K[0][1] = mult * (aN_x*lambda*N_y+aN_y*mu*N_x);
  K[0][2] = mult * (aN_x*lambda*N_z+aN_z*mu*N_x);
  K[1][0] = mult * (aN_y*lambda*N_x+aN_x*mu*N_y);
  K[1][1] = mult * (aN_y*N_y*lambda+2.0*aN_y*mu*N_y+aN_x*mu*N_x+aN_z*mu*N_z);
  K[1][2] = mult * (aN_y*lambda*N_z+aN_z*mu*N_y);
  K[2][0] = mult * (aN_z*lambda*N_x+aN_x*mu*N_z);
  K[2][1] = mult * (aN_z*lambda*N_y+aN_y*mu*N_z);
  K[2][2] = mult * (aN_z*N_z*lambda+2.0*aN_z*mu*N_z+aN_y*mu*N_y+aN_x*mu*N_x);
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
