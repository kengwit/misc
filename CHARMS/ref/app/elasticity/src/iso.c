static void
btcb_3d_aniso_symm_C (double N_x, double N_y, double N_z,
                      double aN_x, double aN_y, double aN_z,
                      double C[6][6], double K[3][3])
{
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
  K[0][0] = N_x*aN_x*C00+N_x*aN_y*C30+N_x*aN_z*C50+N_y*aN_x*C03+N_y*aN_y*C33+N_y*aN_z*C53+N_z*aN_x*C05+N_z*aN_y*C35+N_z*aN_z*C55;
  K[0][1] = N_y*aN_x*C01+N_y*aN_y*C31+N_y*aN_z*C51+N_x*aN_x*C03+N_x*aN_y*C33+N_x*aN_z*C53+N_z*aN_x*C04+N_z*aN_y*C34+N_z*aN_z*C54;
  K[0][2] = N_z*aN_x*C02+N_z*aN_y*C32+N_z*aN_z*C52+N_y*aN_x*C04+N_y*aN_y*C34+N_y*aN_z*C54+N_x*aN_x*C05+N_x*aN_y*C35+N_x*aN_z*C55;
  K[1][0] = N_x*aN_y*C10+N_x*aN_x*C30+N_x*aN_z*C40+N_y*aN_y*C13+N_y*aN_x*C33+N_y*aN_z*C43+N_z*aN_y*C15+N_z*aN_x*C35+N_z*aN_z*C45;
  K[1][1] = N_y*aN_y*C11+N_y*aN_x*C31+N_y*aN_z*C41+N_x*aN_y*C13+N_x*aN_x*C33+N_x*aN_z*C43+N_z*aN_y*C14+N_z*aN_x*C34+N_z*aN_z*C44;
  K[1][2] = N_z*aN_y*C12+N_z*aN_x*C32+N_z*aN_z*C42+N_y*aN_y*C14+N_y*aN_x*C34+N_y*aN_z*C44+N_x*aN_y*C15+N_x*aN_x*C35+N_x*aN_z*C45;
  K[2][0] = N_x*aN_z*C20+N_x*aN_y*C40+N_x*aN_x*C50+N_y*aN_z*C23+N_y*aN_y*C43+N_y*aN_x*C53+N_z*aN_z*C25+N_z*aN_y*C45+N_z*aN_x*C55;
  K[2][1] = N_y*aN_z*C21+N_y*aN_y*C41+N_y*aN_x*C51+N_x*aN_z*C23+N_x*aN_y*C43+N_x*aN_x*C53+N_z*aN_z*C24+N_z*aN_y*C44+N_z*aN_x*C54;
  K[2][2] = N_z*aN_z*C22+N_z*aN_y*C42+N_z*aN_x*C52+N_y*aN_z*C24+N_y*aN_y*C44+N_y*aN_x*C54+N_x*aN_z*C25+N_x*aN_y*C45+N_x*aN_x*C55;
}

static void
btcb_3d_iso_symm_C (double N_x, double N_y, double N_z,
                    double aN_x, double aN_y, double aN_z,
                    double lambda, double mu, double K[3][3])
{
  K[0][0] = (N_x*aN_x*lambda)+(2.0*N_x*mu*aN_x)+(N_y*mu*aN_y)+(N_z*mu*aN_z);
  K[0][1] = K[1][0] = (N_x*lambda*aN_y)+(N_y*mu*aN_x);
  K[0][2] = K[2][0] = (N_x*lambda*aN_z)+(N_z*mu*aN_x);
  K[1][1] = (N_y*aN_y*lambda)+(2.0*N_y*mu*aN_y)+(N_x*mu*aN_x)+(N_z*mu*aN_z);
  K[1][2] = K[2][1] = (N_y*lambda*aN_z)+(N_z*mu*aN_y);
  K[2][2] = (N_z*aN_z*lambda)+(2.0*N_z*mu*aN_z)+(N_y*mu*aN_y)+(N_x*mu*aN_x);
}
