clear all;
close all;
clc;
format short e

% generate random symmetric +ve definite Be (left Cauchy-Green)
%Be=rand(3,3);
%Be = Be+Be';
%Be = Be + 3*eye(3);

% generate random symmetric +ve definite Be (left Cauchy-Green)
Be2D=rand(2,2);
Be2D = Be2D+Be2D';
Be2D = Be2D + 2*eye(2);
Be = zeros(3,3);
Be(1:2,1:2) = [ 1 0
                0 2.2];
Be(3,3) = 1;

Be
[BeVec,BePr]=eig(Be);
BePr=[BePr(1) BePr(5) BePr(9)]';

%
% overwrite with diagonals only
%
%  L = parDerGen(Be,BeVec,BePr,log(BePr),1./BePr);
%
X     = Be;
eV    = BeVec;
eP    = BePr; 
yP    = log(BePr);
ydash = 1./BePr;
L1    = parDerGenWILL_annotated(X,eV,eP,yP,ydash)

L2 = parDerGenKWL(X,eV,eP)

%norm(L1-L2)

L3 = parDerGenKWL2D(X)
