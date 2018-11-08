#include "ParDerGen.h"

/*function L = parDerGenWILL_annotated(X,eV,eP,yP,ydash)
% get eigenvalues of Be (left Cauchy-Green)
%  [BeVec,BePr]=eig(Be);
%
% overwrite with diagonals only
%  BePr=[BePr(1) BePr(5) BePr(9)]'; 
%
%  L = parDerGen(Be,BeVec,BePr,log(BePr),1./BePr);
%
% X = Be
% eV = BeVec
% eP = BePr (diagonals)
% yP = log(BePr)
% ydash = 1./BePr
%*/
void ParDerGen(Eigen::Matrix3d & X, Eigen::MatrixXd & L9)
{
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(X);     // [vec val] = eig(A)
	Eigen::Matrix3d eV = eig.eigenvectors();               // vec
	Eigen::Vector3d eP = eig.eigenvalues();                // diag(val)
	Eigen::Vector3d yP;	for ( int i = 0; i < 3; ++i ) yP(i) = log(eP(i));
	Eigen::Vector3d ydash;	for ( int i = 0; i < 3; ++i ) ydash(i) = 1.0/eP(i);
	
	// just make sure that eigenvalues are not negative
	assert(eP(0) >= 0.0);
	assert(eP(1) >= 0.0);
	assert(eP(2) >= 0.0);
	
	double tol=1e-9; 
	
	// Is=[  eye(3) zeros(3); 
	//	   zeros(3) eye(3)/2]; 
	Eigen::MatrixXd Is;	Is.resize(6,6); Is.fill(0);
	Is.block<3,3>(0,0) = Eigen::MatrixXd::Identity(3,3);    
	Is.block<3,3>(3,3) = 0.5*Is.block<3,3>(0,0);    
	//cout << "Is: " << endl << Is << endl;
	
	// bm1=[1 1 1 0 0 0].';
	Eigen::VectorXd bm1; bm1.resize(6); bm1.fill(0); 
	bm1(0) = 1;
	bm1(1) = 1;
	bm1(2) = 1;
	//cout << "bm1: " << endl << bm1 << endl;
	
	Eigen::MatrixXd L; L.resize(6,6);  
	if (fabs(eP(0))<tol && fabs(eP(1))<tol && fabs(eP(2))<tol) { // all 3 eigenvalues are ZERO
    	cout << "case1\n";
		//all eigenvalues are (close to) zero
		L=Is;
	} else if ( fabs(eP(0)-eP(1))<tol && fabs(eP(0)-eP(2))<tol ) { // all 3 eigenvalues are equal but NOT ZERO
		cout << "case2\n";
		// all eigenvalues are (approx.) the SAME
		// d_ab = d(ln \lam_a)/d(\lam_b) = 1/lam_a * delta_ab
		// and alternative form of outlined by Ogden,
		// 0.5*( ln \lam_a - ln \lam_b )/( \lam_a - \lam_b ) 
		// is replaced by 0.5*d( ln \lam_a - ln \lam_b )/d( \lam_a  ) 
		// i.e. from l'Hospital rule
		// = 0.5*( 1 - delta_ab )/\lam_a = 0.5/\lam_a for a != b	
		// take derivative directly on G matrix
		L=ydash(0)*Is; // note: yd(0)=1/BePr(0);
		
	} else if ( fabs(eP(0)-eP(1))<tol || fabs(eP(1)-eP(2))<tol || fabs(eP(0)-eP(2))<tol ) { // 2 eigenvalues are equal
	      
		double xa,xc,ya,yc,yda,ydc;
		
		if ( fabs(eP(0)-eP(1)) <tol ) {
			xa=eP(2); 
			xc=eP(0); 
			ya=yP(2); 
			yc=yP(0); 
			yda=ydash(2); 
			ydc=ydash(0); 
		} else if ( abs(eP(1)-eP(2))<tol ) {
			xa=eP(0); 
			xc=eP(1); 
			ya=yP(0); 
			yc=yP(1); 
			yda=ydash(0); 
			ydc=ydash(1); 
		} else {                         
			xa=eP(1); 
			xc=eP(0); 
			ya=yP(1); 
			yc=yP(0); 
			yda=ydash(1); 
			ydc=ydash(0); 
		}
		//
		// ya = log(lambda_a)
		// xa = lambda_a
		// dy_a/dx_b = 0 for a not eq. to b
		// dy_a/dx_a = 1/lambda_a for a == b
		// apply this to Eq. A.47 in complas peric
		//
		// yda = dy_a/dx_a
		// ydb = dy_b/dx_b
		// ydc = dy_c/dx_c
		// and cross terms (e.g. dy_a/dx_b) are zero
		//
		// X = 1 4 7 
		//     2 5 8
		//     3 6 9
		Eigen::VectorXd x; x.resize(6);
		// note: X=Be is symmetric
		x << X(0,0), // x(0) = X(1)  
			 X(1,1), // x(1) = X(5)
			 X(2,2), // x(2) = X(9)
			 X(0,1), // x(3) = X(4)=X(2)
			 X(1,2), // x(4) = X(8)=X(6)
			 X(0,2); // x(5) = X(7)=X(3)
			 			 
		Eigen::VectorXd s; s.resize(5); s.fill(0);
			 
		double Dy  = ya-yc;
		double Dx  = xa-xc; 
		double Dx2 = Dx*Dx; 
		double Dx3 = Dx2*Dx; 
		
		s(0)=Dy/Dx2-ydc/Dx;  
		s(1)=2*xc*Dy/Dx2-(xa+xc)/Dx*ydc;
		s(2)=2*Dy/Dx3-(yda+ydc)/Dx2; 
		s(3)=xc*s(2); // eq. to s4 and s5 in A.47
		s(4)=xc*xc*s(2); // eq. to s6 in A.47
		
		Eigen::MatrixXd dX2dX; dX2dX.resize(6,6);
		
		dX2dX << 2*x(0), 0     , 0     , x(3)        , 0            , x(5)         ,
				 0     , 2*x(1), 0     , x(3)        , x(4)         , 0            ,
				 0     , 0     , 2*x(2), 0           , x(4)         , x(5)         ,				 
				 x(3)  , x(3)  , 0     ,(x(0)+x(1))/2, x(5)/2       , x(4)/2       ,
				 0     , x(4)  , x(4)  , x(5)/2      , (x(1)+x(2))/2, x(3)/2       ,
				 x(5)  , 0     , x(5)  , x(4)/2      , x(3)/2       , (x(0)+x(2))/2;
				 
		L=s(0)*dX2dX-s(1)*Is-s(2)*(x*x.transpose())+s(3)*(x*bm1.transpose()+bm1*x.transpose())-s(4)*(bm1*bm1.transpose());

	} else { // distinct eigenvalues
		
		// use's Miehe's method
		// B = Be^tr
		// G(B) = ln(B) = ln(Be^tr)
		// g_a = ln(lambda_a) where lambda_a = eigenvalues of G(B) = ln(Be^tr)
		// d_ab = d(g_a)/d(lambda_b) = 1/lambda_a (if a==b), zero otherwise
		// M^a_ij = 2nd-order tensor(matrix) formed by tensor product of a-th eigenvectors of G(B)
		//         = converted to vector form (eDir) below  
		
		Eigen::Vector3d D; D.resize(3);
		D << (eP(0)-eP(1))*(eP(0)-eP(2)),
		     (eP(1)-eP(0))*(eP(1)-eP(2)),
		     (eP(2)-eP(0))*(eP(2)-eP(1));
		
		Eigen::Vector3d gama; gama.resize(3); gama.fill(0);
		Eigen::MatrixXd eDir; eDir.resize(6,3); eDir.fill(0);
		
		double alfa=0; 
		double bta=0; 
		for ( int i = 0; i < 3; ++i ) 
		{
			alfa=alfa+yP(i)*eP(i)/D(i); 
			bta=bta+yP(i)/D(i)*X.determinant();
			for ( int j = 0; j < 3; ++j ) 
			{
				gama(i)=gama(i)+yP(j)*eP(j)/D(j)*( X.determinant()/eP(j)-eP(i)*eP(i) ) / ( eP(i)*eP(i) );
			}
			
			Eigen::Matrix3d esq = eV.col(i)*eV.col(i).transpose(); 
			// note: esq is symmetric
			eDir(0,i) = esq(0,0); 
			eDir(1,i) = esq(1,1); 
			eDir(2,i) = esq(2,2); 
			eDir(3,i) = esq(0,1); 
			eDir(4,i) = esq(1,2); 
			eDir(5,i) = esq(2,0);		
		}
		
		Eigen::Matrix3d Xinv = X.inverse(); 
		Eigen::VectorXd y; y.resize(6);
		y << Xinv(0,0), // y(0) = Xinv(1)  
			 Xinv(1,1), // y(1) = Xinv(5)
			 Xinv(2,2), // y(2) = Xinv(9)
			 Xinv(0,1), // y(3) = Xinv(4)=Xinv(2)
			 Xinv(1,2), // y(4) = Xinv(8)=Xinv(6)
			 Xinv(0,2); // y(5) = Xinv(7)=Xinv(3)
		
		Eigen::MatrixXd Ib; Ib.resize(6,6); 
		Ib << y(0)*y(0),    y(3)*y(3),    y(5)*y(5),     y(0)*y(3)          ,     y(3)*y(5)          ,     y(0)*y(5)          ,
			  y(3)*y(3),    y(1)*y(1),    y(4)*y(4),     y(1)*y(3)          ,     y(1)*y(4)          ,     y(3)*y(4)          ,
			  y(5)*y(5),    y(4)*y(4),    y(2)*y(2),     y(4)*y(5)          ,     y(2)*y(4)          ,     y(2)*y(5)          ,
			  y(0)*y(3),    y(1)*y(3),    y(4)*y(5), (y(0)*y(1)+y(3)*y(3))/2, (y(3)*y(4)+y(1)*y(5))/2, (y(0)*y(4)+y(3)*y(5))/2,
			  y(3)*y(5),    y(1)*y(4),    y(2)*y(4), (y(3)*y(4)+y(1)*y(5))/2, (y(2)*y(1)+y(4)*y(4))/2, (y(2)*y(3)+y(4)*y(5))/2,
			  y(0)*y(5),    y(3)*y(4),    y(2)*y(5), (y(0)*y(4)+y(3)*y(5))/2, (y(2)*y(3)+y(4)*y(5))/2, (y(2)*y(0)+y(5)*y(5))/2;

		L=alfa*Is-bta*Ib; 
		for ( int i = 0; i < 3; ++i ) 
		{
			L=L+(ydash(i)+gama(i))*eDir.col(i)*eDir.col(i).transpose(); 
		}
			
			
	}
	
	// convert L to a 9x9 version
	L9.resize(9,9);
	L9 << L(0,0), L(0,1), L(0,2), L(0,3), L(0,3), L(0,4), L(0,4), L(0,5), L(0,5),
		  L(1,0), L(1,1), L(1,2), L(1,3), L(1,3), L(1,4), L(1,4), L(1,5), L(1,5),
	      L(2,0), L(2,1), L(2,2), L(2,3), L(2,3), L(2,4), L(2,4), L(2,5), L(2,5),
		  L(3,0), L(3,1), L(3,2), L(3,3), L(3,3), L(3,4), L(3,4), L(3,5), L(3,5),
		  L(3,0), L(3,1), L(3,2), L(3,3), L(3,3), L(3,4), L(3,4), L(3,5), L(3,5),
		  L(4,0), L(4,1), L(4,2), L(4,3), L(4,3), L(4,4), L(4,4), L(4,5), L(4,5),
		  L(4,0), L(4,1), L(4,2), L(4,3), L(4,3), L(4,4), L(4,4), L(4,5), L(4,5),
		  L(5,0), L(5,1), L(5,2), L(5,3), L(5,3), L(5,4), L(5,4), L(5,5), L(5,5),
		  L(5,0), L(5,1), L(5,2), L(5,3), L(5,3), L(5,4), L(5,4), L(5,5), L(5,5);
		  
		  
	
}

