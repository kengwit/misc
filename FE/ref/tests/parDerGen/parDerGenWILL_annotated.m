function L = parDerGenWILL_annotated(X,eV,eP,yP,ydash)
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
%
tol=1e-12; 

Is=[eye(3) zeros(3); 
    zeros(3) eye(3)/2]; 

bm1=[1 1 1 0 0 0].';
if (abs(eP(1))<tol && abs(eP(2))<tol && abs(eP(3))<tol) % all 3 eigenvalues are ZERO
    % all eigenvalues are (close to) zero
	L=Is;
elseif abs(eP(1)-eP(2))<tol && abs(eP(1)-eP(3))<tol % all 3 eigenvalues are equal but NOT ZERO
    
    % all eigenvalues are (approx.) the SAME
    % d_ab = d(ln \lam_a)/d(\lam_b) = 1/lam_a * delta_ab
	% and alternative form of outlined by Ogden,
	% 0.5*( ln \lam_a - ln \lam_b )/( \lam_a - \lam_b ) 
	% is replaced by 0.5*d( ln \lam_a - ln \lam_b )/d( \lam_a  ) 
	% i.e. from l'Hospital rule
	% = 0.5*( 1 - delta_ab )/\lam_a = 0.5/\lam_a for a != b	
	% take derivative directly on G matrix
	L=ydash(1)*Is; % note: yd(1)=1/BePr(1);    
elseif abs(eP(1)-eP(2))<tol || abs(eP(2)-eP(3))<tol || abs(eP(1)-eP(3))<tol % 2 eigenvalues are equal
  if abs(eP(1)-eP(2))<tol
	xa=eP(3); 
	xc=eP(1); 
	ya=yP(3); 
	yc=yP(1); 
	yda=ydash(3); 
	ydc=ydash(1); 
  elseif abs(eP(2)-eP(3))<tol
	xa=eP(1); 
	xc=eP(2); 
	ya=yP(1); 
	yc=yP(2); 
	yda=ydash(1); 
	ydc=ydash(2); 
  else                         
	xa=eP(2); 
	xc=eP(1); 
	ya=yP(2); 
	yc=yP(1); 
	yda=ydash(2); 
	ydc=ydash(1); 
  end 
  % ya = log(lambda_a)
  % xa = lambda_a
  % dy_a/dx_b = 0 for a not eq. to b
  % dy_a/dx_a = 1/lambda_a for a == b
  % apply this to Eq. A.47 in complas peric
  %
  % yda = dy_a/dx_a
  % ydb = dy_b/dx_b
  % ydc = dy_c/dx_c
  % and cross terms (e.g. dy_a/dx_b) are zero
  %
  x=X([1 5 9 4 6 3]).'; % note: X=Be is symmetric
  s=zeros(5,1);
  s(1)=(ya-yc)/(xa-xc)^2-ydc/(xa-xc);  
  s(2)=2*xc*(ya-yc)/(xa-xc)^2-(xa+xc)/(xa-xc)*ydc;
  s(3)=2*(ya-yc)/(xa-xc)^3-(yda+ydc)/(xa-xc)^2; 
  s(4)=xc*s(3); % eq. to s4 and s5 in A.47
  s(5)=xc^2*s(3); % eq. to s6 in A.47
 
  dX2dX=[2*X(1) 0      0      X(2)         0             X(3)          ;
         0      2*X(5) 0      X(2)         X(6)          0             ;
         0      0      2*X(9) 0            X(6)          X(3)          ;
         X(2)   X(2)   0     (X(1)+X(5))/2 X(3)/2        X(6)/2        ;
         0      X(6)   X(6)   X(3)/2       (X(5)+X(9))/2 X(2)/2        ;
         X(3)   0      X(3)   X(6)/2       X(2)/2        (X(1)+X(9))/2];
  L=s(1)*dX2dX-s(2)*Is-s(3)*(x*x.')+s(4)*(x*bm1.'+bm1*x.')-s(5)*(bm1*bm1.');
else % distinct eigenvalues
  % use's Miehe's method
  % B = Be^tr
  % G(B) = ln(B) = ln(Be^tr)
  % g_a = ln(lambda_a) where lambda_a = eigenvalues of G(B) = ln(Be^tr)
  % d_ab = d(g_a)/d(lambda_b) = 1/lambda_a (if a==b), zero otherwise
  % M^a_ij = 2nd-order tensor(matrix) formed by tensor product of a-th eigenvectors of G(B)
  %         = converted to vector form (eDir) below  
  D=[(eP(1)-eP(2))*(eP(1)-eP(3));
     (eP(2)-eP(1))*(eP(2)-eP(3));
	 (eP(3)-eP(1))*(eP(3)-eP(2))];
  
  alfa=0; 
  bta=0; 
  gama=zeros(3,1); 
  eDir=zeros(6,3);
  
  for i=1:3 
    alfa=alfa+yP(i)*eP(i)/D(i); 
	bta=bta+yP(i)/D(i)*det(X);
    for j=1:3
		gama(i)=gama(i)+yP(j)*eP(j)/D(j)*(det(X)/eP(j)-eP(i)^2)*1/eP(i)^2;
	end
    esq=eV(:,i)*eV(:,i).'; 
	eDir(:,i)=[esq(1,1) esq(2,2) esq(3,3) esq(1,2) esq(2,3) esq(3,1)].';
  end
  y=inv(X);
  Ib=[y(1)^2    y(2)^2    y(7)^2     y(1)*y(2)               y(2)*y(7)               y(1)*y(7)              ;
      y(2)^2    y(5)^2    y(6)^2     y(5)*y(2)               y(5)*y(6)               y(2)*y(6)              ;
      y(7)^2    y(6)^2    y(9)^2     y(6)*y(7)               y(9)*y(6)               y(9)*y(7)              ;
      y(1)*y(2) y(5)*y(2) y(6)*y(7) (y(1)*y(5)+y(2)^2)/2    (y(2)*y(6)+y(5)*y(7))/2 (y(1)*y(6)+y(2)*y(7))/2 ;
      y(2)*y(7) y(5)*y(6) y(9)*y(6) (y(2)*y(6)+y(5)*y(7))/2 (y(9)*y(5)+y(6)^2)/2    (y(9)*y(2)+y(6)*y(7))/2 ;
      y(1)*y(7) y(2)*y(6) y(9)*y(7) (y(1)*y(6)+y(2)*y(7))/2 (y(9)*y(2)+y(6)*y(7))/2 (y(9)*y(1)+y(7)^2)/2   ];
  
  L=alfa*Is-bta*Ib; 
  for i=1:3
	L=L+(ydash(i)+gama(i))*eDir(:,i)*eDir(:,i).'; 
  end
end

% convert L to a 9x9 version
L=[          L(1:3,1:3)           L(1:3,[4 4 5 5 6 6]);
   L([4 4 5 5 6 6],1:3) L([4 4 5 5 6 6],[4 4 5 5 6 6])];