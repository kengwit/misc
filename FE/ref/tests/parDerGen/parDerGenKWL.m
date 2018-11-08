function L = parDerGenKWL(X,eV,eP)
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
    yP = log(eP);
    ydash = 1./eP;
    
    tol=1.e-12; 
    
    L = zeros(3,3,3,3);
    
    Iden=eye(3);    

    Is = zeros(3,3,3,3);
    for i=1:3
        for j=1:3
            for k=1:3
                for l=1:3
                   Is(i,j,k,l) = 0.5*(Iden(i,k)*Iden(j,l)+Iden(i,l)*Iden(j,k));
                end
            end     
        end
    end


    if (abs(eP(1))<tol && abs(eP(2))<tol && abs(eP(3))<tol) % all 3 eigenvalues are ZERO
        fprintf('3 eigenvalues are close to zero\n');
        % all eigenvalues are (close to) zero
        L=Is;
    elseif abs(eP(1)-eP(2))<tol && abs(eP(1)-eP(3))<tol % all 3 eigenvalues are equal but NOT ZERO
        fprintf('3 eigenvalues are equal but not zero\n');
        
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
        fprintf('2 eigenvalues are equal\n');
        
        if ( abs(eP(1)-eP(2))<tol )
            % a=3 != b=1,c=2
            xa=eP(3); 
            xc=eP(1); 
            ya=yP(3); 
            yc=yP(1); 
            yda=ydash(3); 
            ydc=ydash(1); 
        elseif ( abs(eP(2)-eP(3))<tol )
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
    %   % ya = log(lambda_a)
    %   % xa = lambda_a
    %   % dy_a/dx_b = 0 for a not eq. to b
    %   % dy_a/dx_a = 1/lambda_a for a == b
    %   % apply this to Eq. A.47 in complas peric
    %   %
    %   % yda = dy_a/dx_a
    %   % ydb = dy_b/dx_b
    %   % ydc = dy_c/dx_c
    %   % and cross terms (e.g. dy_a/dx_b) are zero
    %   %

        dX2dX = zeros(3,3,3,3);
        XxX   = zeros(3,3,3,3);
        XxI   = zeros(3,3,3,3);
        IxX   = zeros(3,3,3,3);
        IxI   = zeros(3,3,3,3);

        for i=1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                       dX2dX(i,j,k,l) = 0.5*(Iden(i,k)*X(l,j)+Iden(i,l)*X(k,j)+Iden(j,l)*X(i,k)+Iden(k,j)*X(i,l));
                       XxX(i,j,k,l) = X(i,j)*X(k,l);
                       XxI(i,j,k,l) = X(i,j)*Iden(k,l);
                       IxX(i,j,k,l) = Iden(i,j)*X(k,l);
                       IxI(i,j,k,l) = Iden(i,j)*Iden(k,l);                  
                    end
                end     
            end
        end


        % see Eq. A.53
        s=zeros(6,1);
        s(1)=(ya-yc)/(xa-xc)^2-ydc/(xa-xc);  
        s(2)=2*xc*(ya-yc)/(xa-xc)^2-(xa+xc)/(xa-xc)*ydc;
        s(3)=2*(ya-yc)/(xa-xc)^3-(yda+ydc)/(xa-xc)^2; 
        s(4)=xc*s(3);   % eq. to s4 and s5 in A.47
        s(5)=s(4);      % eq. to s4 and s5 in A.47
        s(6)=xc^2*s(3); % eq. to s6 in A.47

        L = s(1)*dX2dX - s(2)*Is - s(3)*XxX + s(4)*XxI + s(5)*IxX - s(6)*IxI;

    else % distinct eigenvalues
        fprintf('distinct eigenvalues\n');
        dX2dX = zeros(3,3,3,3);
        for i=1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                       dX2dX(i,j,k,l) = 0.5*(Iden(i,k)*X(l,j)+Iden(i,l)*X(k,j)+Iden(j,l)*X(i,k)+Iden(k,j)*X(i,l));
                    end
                end     
            end
        end

        % eigenprojections
        E=zeros(3,3,3);
        E(1,:,:) = eV(:,1)*eV(:,1)';
        E(2,:,:) = eV(:,2)*eV(:,2)';
        E(3,:,:) = eV(:,3)*eV(:,3)';
        
        
        EaxEa = zeros(3,3,3,3);
        EbxEb = zeros(3,3,3,3);
        EcxEc = zeros(3,3,3,3);
        
        permut=[1 2 3;
                2 3 1;
                3 1 2];
        
        for p=1:3
            
            id=permut(p,:);
            i1 = id(1);
            i2 = id(2);
            i3 = id(3);
            
            for i=1:3
                for j=1:3
                    for k=1:3
                        for l=1:3
                           EaxEa(i,j,k,l) = E(i1,i,j)*E(i1,k,l);
                           EbxEb(i,j,k,l) = E(i2,i,j)*E(i2,k,l);
                           EcxEc(i,j,k,l) = E(i3,i,j)*E(i3,k,l);
                        end
                    end     
                end
            end

            xa=eP(i1); 
            xb=eP(i2); 
            xc=eP(i3); 
            ya=yP(i1); 
            yda=ydash(i1); 

            L=L+ya*( ...
                    dX2dX ...
                  - (xb+xc)*Is ...
                  - ((xa-xb)+(xa-xc))*EaxEa ...
                  - (xb-xc)*(EbxEb - EcxEc) ...
                   )/((xa-xb)*(xa-xc)) + yda*EaxEa;

        end 
        
    end

    % convert L to a 9x9 version
    L=convertToMat9(L);

end

