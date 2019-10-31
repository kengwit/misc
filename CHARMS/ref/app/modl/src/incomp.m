% Pressure/volumetric strain interpolation matrix
Nv=[sym('1-r-s-t') sym('r') sym('s') sym('t')];
% Integration point locations
gpcoord = [ 0.13819660 0.13819660 0.13819660 0.58541020 ;
            0.58541020 0.13819660 0.13819660 0.13819660 ;
            0.13819660 0.58541020 0.13819660 0.13819660 ;
            0.13819660 0.13819660 0.58541020 0.13819660 ];
% Quadrature weight
gpw = 0.041666666666666666667;

% E matrix
E=zeros(4,4);
for gp = 1:4
   r = gpcoord(gp,1);
   s = gpcoord(gp,2);
   t = gpcoord(gp,3);
   u = gpcoord(gp,4);
   Nveval = eval (Nv);
   E = E + Nveval' * Nveval * gpw;
end

% B matrix at the barycenter
r=0.25; s=0.25; t=0.25; u=0.25;
B1 = [ diff(Nv(1),'r') 0               0 ; 
       0               diff(Nv(1),'s') 0;
       0               0               diff(Nv(1),'t');
       diff(Nv(1),'s') diff(Nv(1),'r') 0;
       0               diff(Nv(1),'t') diff(Nv(1),'s');
       diff(Nv(1),'t') 0               diff(Nv(1),'r') ];
B2 = [ diff(Nv(2),'r') 0               0 ; 
       0               diff(Nv(2),'s') 0;
       0               0               diff(Nv(2),'t');
       diff(Nv(2),'s') diff(Nv(2),'r') 0;
       0               diff(Nv(2),'t') diff(Nv(2),'s');
       diff(Nv(2),'t') 0               diff(Nv(2),'r') ];
B3 = [ diff(Nv(3),'r') 0               0 ; 
       0               diff(Nv(3),'s') 0;
       0               0               diff(Nv(3),'t');
       diff(Nv(3),'s') diff(Nv(3),'r') 0;
       0               diff(Nv(3),'t') diff(Nv(3),'s');
       diff(Nv(3),'t') 0               diff(Nv(3),'r') ];
B4 = [ diff(Nv(4),'r') 0               0 ; 
       0               diff(Nv(4),'s') 0;
       0               0               diff(Nv(4),'t');
       diff(Nv(4),'s') diff(Nv(4),'r') 0;
       0               diff(Nv(4),'t') diff(Nv(4),'s');
       diff(Nv(4),'t') 0               diff(Nv(4),'r') ];
B = [ B1 B2 B3 B4 ] ;

m=[1 1 1 0 0 0]';
C=B'*m*eval(Nv);

W=eval(inv(E)*C');
Id = (eye(6)-1/3*m*m');
I0 = 1/2*[ 2 0 0 0 0 0 ; 
           0 2 0 0 0 0 ;
           0 0 2 0 0 0 ;
           0 0 0 1 0 0 ;
           0 0 0 0 1 0 ;
           0 0 0 0 0 1 ];
Bbar = Id * B + 1/3 * m * eval(Nv)*W;
