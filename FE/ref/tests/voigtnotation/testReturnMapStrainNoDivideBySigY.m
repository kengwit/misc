clear all;
close all;
clc;
format long
% =======================================
tol = 1.e-12;
% 9x6 Projection matrix
Pmat = [1 0 0   0   0   0;
        0 1 0   0   0   0;
        0 0 1   0   0   0;
        0 0 0 0.5   0   0;
        0 0 0 0.5   0   0;
        0 0 0   0 0.5   0;
        0 0 0   0 0.5   0;
        0 0 0   0   0 0.5;
        0 0 0   0   0 0.5];
    
bm1 = zeros(9,1);
bm1(1:3)=1.0;

Idev = eye(9) - 1./3.*bm1*bm1';
% =======================================
% cauchy stress tensor
E = 1.0;
nu = 0.3;
sigmay=0.01;
C11 = E*(1-nu)/((1+nu)*(1-2*nu));
C12 = E*nu/((1+nu)*(1-2*nu));
G2 = C11-C12;
De9 = [C11 C12 C12  0  0  0  0  0  0;
      C12 C11 C12  0  0  0  0  0  0;
      C12 C12 C11  0  0  0  0  0  0;
        0   0   0 G2  0  0  0  0  0;
        0   0   0  0 G2  0  0  0  0;  
        0   0   0  0  0 G2  0  0  0;  
        0   0   0  0  0  0 G2  0  0;    
        0   0   0  0  0  0  0 G2  0;  
        0   0   0  0  0  0  0  0 G2];

Gmod = E/(2.0*(1.0+nu));     % shear modulus
Kmod = E/(3.0*(1.0-2.0*nu)); % bulk modulus
C11 = Kmod+4.0*Gmod/3.0;
C12 = Kmod-2.0*Gmod/3.0;
G2  = C11-C12;
De2 = [C11 C12 C12  0  0  0  0  0  0;
      C12 C11 C12  0  0  0  0  0  0;
      C12 C12 C11  0  0  0  0  0  0;
        0   0   0 G2  0  0  0  0  0;
        0   0   0  0 G2  0  0  0  0;  
        0   0   0  0  0 G2  0  0  0;  
        0   0   0  0  0  0 G2  0  0;    
        0   0   0  0  0  0  0 G2  0;  
        0   0   0  0  0  0  0  0 G2]; 
Ce9=inv(De9);
% =======================================

% sigxx = rand;
% sigyy = rand;
% sigxy = rand;
% sigzz = nu*(sigxx+sigyy);
% sig = [sigxx sigxy     0;
%        sigxy sigyy     0;
%            0     0 sigzz];

epsEtrxx = 0.1;
epsEtryy = 0.1;
epsEtrxy = 0.05;
epsEtr = [epsEtrxx epsEtrxy   0;
          epsEtrxy epsEtryy   0;
                 0        0   0];
epsEtr9    = zeros(9,1);
epsEtr9(1) = epsEtr(1,1);
epsEtr9(2) = epsEtr(2,2);
epsEtr9(3) = epsEtr(3,3);
epsEtr9(4) = epsEtr(1,2);
epsEtr9(5) = epsEtr(2,1);
epsEtr9(6) = epsEtr(2,3);
epsEtr9(7) = epsEtr(3,2);
epsEtr9(8) = epsEtr(3,1);
epsEtr9(9) = epsEtr(1,3);

% set
epsE9 = epsEtr9;

% trial sigma
Sig9 = De9*epsEtr9;

% deviatoric stress tensor
trsig = (Sig9(1)+Sig9(2)+Sig9(3));

sdev9 = Sig9 - trsig/3.0*bm1;
sdevnorm = norm(sdev9);

j2 = 0.5*sdevnorm*sdevnorm;
	
% yield function value
f = sqrt(3.0*j2)-sigmay;
f_Borja = sdevnorm-sqrt(2.0/3.0)*sigmay;
sdev9_Borja = sdev9;
if ( f > 0.0 )
    fprintf('yield\n');
    
    maxit = 5;
	dgam = 0.0;
	converged_flag = false;
	
    % initialize residual vector
	resid = zeros(10,1);
	resid(10) = f;
			
	nbar = sdev9/sdevnorm;	
	df  = sqrt(3.0/2.0)*nbar;
	ddf = sqrt(3)/(2.0*sqrt(j2))*( Idev - nbar*nbar' ); 
   
    % newton iteration
	for itnum = 1:maxit
		
		fprintf('\tnewton iter = %d\n',itnum);
            
        % form Jacobian
			
			
        A = zeros(10,10);
        A(1:9,1:9) = eye(9)+dgam*ddf*De9;
        A(10,1:9)  = df'*De9;
        A(1:9,10)  = df;

        
        
        Ainv = inv(A);
			
		dx = Ainv * resid;
		
		epsE9 = epsE9 - dx(1:9); 			
		dgam  = dgam  - dx(10);      
		
        Sig9 = De9*epsE9
		
        % deviatoric stress tensor
        trsig = (Sig9(1)+Sig9(2)+Sig9(3));

        sdev9 = Sig9 - trsig/3.0*bm1;
        sdevnorm = norm(sdev9);

        j2 = 0.5*sdevnorm*sdevnorm
        
		% recompute terms in jacobian matrix for next iteration
		nbar = sdev9/sdevnorm;	
        df  = sqrt(3.0/2.0)*nbar
        
        %nbar*nbar'
        %sqrt(3)/(2.0*sqrt(j2))
        %Idev
        
        ddf = sqrt(3)/(2.0*sqrt(j2))*( Idev - nbar*nbar' ) 
   
        % recompute residual vector
		resid = [ epsE9 - epsEtr9 + dgam*df;
                  sqrt(3.0*j2)-sigmay ];
              
        fprintf('norm(resid(1:9)) = %e, abs(resid(10)) = %e \n',norm(resid(1:9)),abs(resid(10)));
        if ( ( norm(resid(1:9)) < tol ) && ( abs(resid(10)) < tol ) )
            converged_flag = 1;
            break;
        end

			
    end
		
    if ( converged_flag == true ) 

        % consistent tangent operator
        df
        ddf
        dgam
        A = zeros(10,10);
        A(1:9,1:9) = Ce9+dgam*ddf;
        A(10,1:9)  = df';
        A(1:9,10)  = df;
        
        B = inv(A);
        
        Dalg9 = B(1:9,1:9); 
        Dalg6 = Pmat'*Dalg9*Pmat;
        
        fprintf('SigXX = %f\n',Sig9(1));
        fprintf('SigYY = %f\n',Sig9(2));
        fprintf('SigZZ = %f\n',Sig9(3));
        fprintf('SigXY = %f\n',Sig9(4));
        fprintf('SigYZ = %f\n',Sig9(6));
        fprintf('SigZX = %f\n',Sig9(9));
        
        Dalg6
        
        % Borja style
        sdevnorm_Borja = norm(sdev9_Borja);
        nbar_Borja = sdev9_Borja/sdevnorm_Borja;
        A0 = sdevnorm_Borja;
        
        bm6 = zeros(6,1);
        bm6(1:3)=1.0;

        IDMat0=eye(6);
        IDMat0(4,4)=IDMat0(4,4)/2;
        IDMat0(5,5)=IDMat0(5,5)/2;
        IDMat0(6,6)=IDMat0(6,6)/2;
        Pdev=IDMat0-bm6*bm6'/3;  
        H_Borja=0;
        NTemp = [nbar_Borja(1);
                 nbar_Borja(2);
                 nbar_Borja(3);
                 nbar_Borja(4);
                 nbar_Borja(6);
                 nbar_Borja(9)];
        bulk=Kmod;
        mu=Gmod;
        dlam= f_Borja/(2*mu+H_Borja);
        % Consistent Tangent Operator (Dalg)
        cep   = bulk*(bm6*bm6') + 2*mu*Pdev - ...
               (2*mu/(1+H_Borja/(2*mu)))*(NTemp*NTemp');
        gamma = 2*mu*dlam/A0;    
        DalgBorja  = cep - 2*mu*gamma*( Pdev - (NTemp*NTemp') )

        norm(Dalg6-DalgBorja)
    else 
        cout << "no convergence in return map after " << maxit << " iterations " << endl;

    end

else
    fprintf('elastic update\n');
end



        

