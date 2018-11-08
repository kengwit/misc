function L = parDerGenKWL2D(X)
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
    tol=1.e-12;
    
    Iden=eye(3);
    Iden(3,3)=0;

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


    % in-plane
    Xinplane=X(1:2,1:2);
    [eV,eP]=eig(Xinplane);    
    
    xP=zeros(3,1);
    yP=zeros(3,1);
    ydash=zeros(3,1);
    
    xP = [eP(1,1); eP(2,2); 1];
    yP = log(xP);
    ydash = 1./xP;
    
    
    % eigenprojections
    E=zeros(2,3,3);
    E(1,1:2,1:2) = eV(:,1)*eV(:,1)';
    E(2,1:2,1:2) = eV(:,2)*eV(:,2)';
    
    eV3=[0;0;1];
    E(3,:,:) = eV3*eV3';
    
    L = zeros(3,3,3,3);
    
    if ( abs(xP(1)-xP(2))<tol ) 
        fprintf('2d same eigevalue\n');
        L=ydash(1)*Is; 
        
    else
        fprintf('2d distinct eigevalue\n');
        
        E1xE1=zeros(3,3,3,3);
        E2xE2=zeros(3,3,3,3);
        
        for i=1:3
            for j=1:3
                for k=1:3
                    for l=1:3
                       E1xE1(i,j,k,l) = E(1,i,j)*E(1,k,l);
                       E2xE2(i,j,k,l) = E(2,i,j)*E(2,k,l);
                       E3xE3(i,j,k,l) = E(3,i,j)*E(3,k,l);
                    end
                end     
            end
        end

        L=(yP(1)-yP(2))*( Is - E1xE1 - E2xE2 ) / (xP(1)-xP(2))...
             + ydash(1)*E1xE1+ydash(2)*E2xE2;
    
    end
    
    % add out-of-plane terms
    L=L+ydash(3)*E3xE3;

    % convert L to a 9x9 version
    L=convertToMat9(L);

end

