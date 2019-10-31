% Compute the B-bar matrix for a tetrahedral element
function [Bbar_mult] = bbar 
% Pressure/volumetric strain interpolation matrix
% Integration point locations
gpcoord = [ 0.13819660 0.13819660 0.13819660 0.58541020 ;
            0.58541020 0.13819660 0.13819660 0.13819660 ;
            0.13819660 0.58541020 0.13819660 0.13819660 ;
            0.13819660 0.13819660 0.58541020 0.13819660 ];
% Quadrature weight
gpw = 0.041666666666666666667;

% Ehat matrix (int(N'*N) over the reference volume)
Ehat=zeros(4,4);
for gp = 1:4
   r = gpcoord(gp,1);
   s = gpcoord(gp,2);
   t = gpcoord(gp,3);
   u = gpcoord(gp,4);
   Nv = Nvp (r, s, t);
   Ehat = Ehat + Nv' * Nv * gpw;
end
Ehat
inv(Ehat)

% B matrix at the barycenter
r=0.25; s=0.25; t=0.25; u=0.25;
Nvp_v = Nvp (r, s, t);

m=[1 1 1 0 0 0]';
% C=B'*m*Nvp_v;

Id = (eye(6)-1/3*m*m');
I0 = 1/2*[ 2 0 0 0 0 0 ; 
           0 2 0 0 0 0 ;
           0 0 2 0 0 0 ;
           0 0 0 1 0 0 ;
           0 0 0 0 1 0 ;
           0 0 0 0 0 1 ];
Bbar_mult = (Id + 1/3 * m * Nvp_v * inv(Ehat) * Nvp_v' * m');
return;

function [symbN] = symb_Nvp 
symbN = [sym('1-r-s-t') sym('r') sym('s') sym('t')];
return

function [evalN] = Nvp (r, s, t)
Nv=symb_Nvp;
evalN = eval(Nv);
return

function [evalNr] = diff_Nvp_r (r, s, t)
Nv=symb_Nvp;
evalNr = eval(diff(Nv,'r'));
return

function [evalNs] = diff_Nvp_s (r, s, t)
Nv=symb_Nvp;
evalNs = eval(diff(Nv,'s'));
return

function [evalNt] = diff_Nvp_t (r, s, t)
Nv=symb_Nvp;
evalNt = eval(diff(Nv,'t'));
return
