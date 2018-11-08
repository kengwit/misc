clear all;
close all;
clc;

% cauchy stress tensor
sig = rand(3);
sig = 0.5*(sig+sig');

% deviatoric stress tensor
trsig = sig(1,1)+sig(2,2)+sig(3,3);
sdev = sig - trsig/3.0*eye(3);
sdevnorm = norm(sdev,'fro');

% flow vector
nbar = sdev/sdevnorm;

% 4th order nbar_ij nbar_kl
nxn = zeros(3,3,3,3);

for i=1:3
for j=1:3
for k=1:3
for l=1:3
    nxn(i,j,k,l) = nbar(i,j)*nbar(k,l);
end
end
end
end

% reduce to 6x6
nxn6 = zeros(6,6);
%
nxn6(1,1) = nxn(1,1,1,1); 
nxn6(1,2) = nxn(1,1,2,2);
nxn6(1,3) = nxn(1,1,3,3);
nxn6(1,4) = nxn(1,1,1,2);
nxn6(1,5) = nxn(1,1,2,3);
nxn6(1,6) = nxn(1,1,1,3);
%
nxn6(2,1) = nxn(2,2,1,1); 
nxn6(2,2) = nxn(2,2,2,2);
nxn6(2,3) = nxn(2,2,3,3);
nxn6(2,4) = nxn(2,2,1,2);
nxn6(2,5) = nxn(2,2,2,3);
nxn6(2,6) = nxn(2,2,1,3);
%
nxn6(3,1) = nxn(3,3,1,1); 
nxn6(3,2) = nxn(3,3,2,2);
nxn6(3,3) = nxn(3,3,3,3);
nxn6(3,4) = nxn(3,3,1,2);
nxn6(3,5) = nxn(3,3,2,3);
nxn6(3,6) = nxn(3,3,1,3);
%
nxn6(4,1) = nxn(1,2,1,1); 
nxn6(4,2) = nxn(1,2,2,2);
nxn6(4,3) = nxn(1,2,3,3);
nxn6(4,4) = nxn(1,2,1,2);
nxn6(4,5) = nxn(1,2,2,3);
nxn6(4,6) = nxn(1,2,1,3);
%
nxn6(5,1) = nxn(2,3,1,1); 
nxn6(5,2) = nxn(2,3,2,2);
nxn6(5,3) = nxn(2,3,3,3);
nxn6(5,4) = nxn(2,3,1,2);
nxn6(5,5) = nxn(2,3,2,3);
nxn6(5,6) = nxn(2,3,1,3);
%
nxn6(6,1) = nxn(1,3,1,1); 
nxn6(6,2) = nxn(1,3,2,2);
nxn6(6,3) = nxn(1,3,3,3);
nxn6(6,4) = nxn(1,3,1,2);
nxn6(6,5) = nxn(1,3,2,3);
nxn6(6,6) = nxn(1,3,1,3);
nxn6
% size 9 approach
sdev9    = zeros(9,1);
sdev9(1) = sdev(1,1);
sdev9(2) = sdev(2,2);
sdev9(3) = sdev(3,3);
sdev9(4) = sdev(1,2);
sdev9(5) = sdev(2,1);
sdev9(6) = sdev(2,3);
sdev9(7) = sdev(3,2);
sdev9(8) = sdev(3,1);
sdev9(9) = sdev(1,3);

nbar9 = sdev9/norm(sdev9);
%fprintf('norm(nbar9) = %f\n',norm(nbar9));

nxn_9x9 = nbar9*nbar9';

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
        

nxn6b = Pmat'*nxn_9x9*Pmat

% coombs approach is incorrect as shown
sdev6    = zeros(6,1);
sdev6(1) = sdev(1,1);
sdev6(2) = sdev(2,2);
sdev6(3) = sdev(3,3);
sdev6(4) = 2*sdev(1,2);
sdev6(5) = 2*sdev(2,3);
sdev6(6) = 2*sdev(1,3);

sdev6 = sdev6/sdevnorm;

nxn6c = sdev6*sdev6'

