function [Sol, NSnormIndicator] = hyperelasticMaterial(solver, solid, pmc, Sol, cnn, crd, ...
                                               elemType, ndof, nen, nElem, BCCyl) %BCStructure BCBarFix

% Quadrature rules for elements
if strcmp(elemType,'3Tri')
    gP = ...
   [1/3,  1/3
    4/3,  1/3
    1/3,  4/3] ;
    gW = ...
   [2/3,  2/3,  2/3] ;
 
    N(:,1) = 0.5.*(2.-gP(:,1)-gP(:,2)) ;
    N(:,2) = 0.5.*(gP(:,1)) ;
    N(:,3) = 0.5.*(gP(:,2)) ;
    
    Nx(:,1) = -0.5.*ones(3,1) ;
    Nx(:,2) =  0.5.*ones(3,1) ;
    Nx(:,3) =  zeros(3,1) ; 
    Ny(:,1) = -0.5.*ones(3,1) ;
    Ny(:,2) =  zeros(3,1) ;
    Ny(:,3) =  0.5.*ones(3,1) ;    
elseif strcmp(elemType,'4Quad')
    gP = ...
   [-5.7735026918962584E-01, -5.7735026918962584E-01
     5.7735026918962584E-01, -5.7735026918962584E-01
    -5.7735026918962584E-01,  5.7735026918962584E-01
     5.7735026918962584E-01,  5.7735026918962584E-01] ;
    gW = [1, 1, 1, 1 ] ;
    
    N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;
    N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
    N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
    N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
    
    Nx(:,1) = -0.25.*(1-gP(:,2)) ;
    Nx(:,2) =  0.25.*(1-gP(:,2)) ;
    Nx(:,3) =  0.25.*(1+gP(:,2)) ;
    Nx(:,4) = -0.25.*(1+gP(:,2)) ;
    Ny(:,1) = -0.25.*(1-gP(:,1)) ;
    Ny(:,2) = -0.25.*(1+gP(:,1)) ;
    Ny(:,3) =  0.25.*(1+gP(:,1)) ;
    Ny(:,4) =  0.25.*(1-gP(:,1)) ;
end

Nx = Nx' ;
Ny = Ny' ;
nQuad = length(gW) ;

velS = velSS;
velSPrev = velSSPrev;
velSDot = velSSDot;
velSDotPrev = velSSDotPrev;
dispS = Sol.dispS;
dispSprev = Sol.dispSPrev;

iif = zeros(nen^2*nElem,1); 
jjf = zeros(nen^2*nElem,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      iif(index+1:index+nElem) = double(cnn(:,i)); 
      jjf(index+1:index+nElem) = double(cnn(:,j));  
      index = index + nElem;
   end
end
        
% Satisfy boundary conditions
bcFix = unique(BCBarFix(:));
velS(bcFix,1,1) = solid.DirichletUval;
velS(bcFix,2,1) = solid.DirichletVval;
 
 
% Interpolate for alpha values for Gen-alpha
velSAlpha = velSPrev + pmc.alpha.*(velS - velSPrev) ;
velSDotAlpha = velSDotPrev + pmc.alphaM.*(velSDot - velSDotPrev) ;
dispSAlpha = dispSPrev + pmc.alpha.*(dispSVel - dispSPrev) ;
        
% Structure equations
xxf = zeros(size(cnn));
yyf = zeros(size(cnn));
ux = zeros(size(cnn));
uy = zeros(size(cnn));
dispx = zeros(size(cnn));
dispy = zeros(size(cnn));

for i=1:nen
   xxf(:,i) = crd(cnn(:,i),1);
   yyf(:,i) = crd(cnn(:,i),2);
   ux(:,i) =  velSAlpha(cnn(:,i),1,1) ;
   uxDot(:,i) = velSDotAlpha(cnn(:,i),1,1) ;
   uy(:,i) =  velSAlpha(cnn(:,i),2,1) ;
   uyDot(:,i) = velSDotAlpha(cnn(:,i),2,1) ;
   dispx(:,i) = dispSAlpha(cnn(:,i),1);
   dispy(:,i) = dispSAlpha(cnn(:,i),2);
   dispxPrev(:,i) = dispSPrev(cnn(:,i),1);
   dispyPrev(:,i) = dispSPrev(cnn(:,i),2);
end
        
% Form element matrix and assemble Galerkin terms

sA = [] ;
sA1 = zeros(nen^2*nElem,nQuad); 
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);

for p = 1:nQuad  
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
           
    volume = abs(volume);
    
    negJacobian = find(volume<0);
    if ~isempty(negJacobian)
       disp('Mesh deformed, Negative Jacobian');
       exit
    end

    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
    locUX  = sum(repmat(N(p,:),nElem,1).*(ux),2);
    locUY  = sum(repmat(N(p,:),nElem,1).*(uy),2);   
    locUXDot  = sum(repmat(N(p,:),nElem,1).*(uxDot),2);
    locUYDot  = sum(repmat(N(p,:),nElem,1).*(uyDot),2);  
    locdispX  = sum(repmat(N(p,:),nElem,1).*(dispx),2);
    locdispY  = sum(repmat(N(p,:),nElem,1).*(dispy),2); 

    locgradXdispX = sum(DNDx.*dispx,2);
    locgradYdispX = sum(DNDy.*dispx,2);
    locgradXdispY = sum(DNDx.*dispy,2);
    locgradYdispY = sum(DNDy.*dispy,2);

    locgradXdispXPrev = sum(DNDx.*dispxPrev,2);
    locgradYdispXPrev = sum(DNDy.*dispxPrev,2);
    locgradXdispYPrev = sum(DNDx.*dispyPrev,2);
    locgradYdispYPrev = sum(DNDy.*dispyPrev,2);

    nu_s = solid.nu;
    % miu_s = YME/(2*(1+niu_s));
    lambda_s = solid.lambda
    % (niu_s*YME)/((1+niu_s)*(1-2niu_s))

    F11  = 1 + locgradXdispXPrev;
    F12  = locgradYdispXPrev;
    F21  = locgradXdispYPrev;
    F22  = 1 + locgradYdispYPrev;

    F11_new  = 1 + locgradXdispX;
    F12_new  = locgradYdispX;
    F21_new  = locgradXdispY;
    F22_new  = 1 + locgradYdispY;

    E11  = 0.5 .* (locgradXdispXPrev + locgradXdispXPrev + locgradXdispXPrev .* locgradXdispXPrev + locgradXdispYPrev .* locgradXdispYPrev);
    E12  = 0.5 .* (locgradYdispXPrev + locgradXdispYPrev + locgradXdispXPrev .* locgradXdispYPrev + locgradXdispYPrev .* locgradXdispYPrev);
    E21  = 0.5 .* (locgradXdispYPrev + locgradYdispXPrev + locgradYdispXPrev .* locgradXdispXPrev + locgradXdispYPrev .* locgradXdispYPrev);
    E22  = 0.5 .* (locgradYdispYPrev + locgradYdispYPrev + locgradYdispXPrev .* locgradXdispYPrev + locgradXdispYPrev .* locgradXdispYPrev);

    E11_new  = 0.5 .* (locgradXdispX + locgradXdispX + locgradXdispX .* locgradXdispX + locgradXdispY .* locgradXdispY);
    E12_new  = 0.5 .* (locgradYdispX + locgradXdispY + locgradXdispX .* locgradXdispY + locgradXdispY .* locgradXdispY);
    E21_new  = 0.5 .* (locgradXdispY + locgradYdispX + locgradYdispX .* locgradXdispX + locgradXdispY .* locgradXdispY);
    E22_new  = 0.5 .* (locgradYdispY + locgradYdispY + locgradYdispX .* locgradXdispY + locgradXdispY .* locgradXdispY);

    TraE = 0.5 .* (2 .* locgradXdispXPrev + 2.* locgradYdispYPrev + ...
                        locgradXdispXPrev    .* locgradXdispXPrev + ...
                        locgradXdispYPrev    .* locgradXdispYPrev + ...
                        locgradYdispXPrev    .* locgradYdispXPrev + ...
                        locgradYdispYPrev    .* locgradYdispYPrev);

    TraE_new = 0.5 .* (2 .* locgradXdispX + 2.* locgradYdispY + ...
                            locgradXdispX    .* locgradXdispX + ...
                            locgradXdispY    .* locgradXdispY + ...
                            locgradYdispX    .* locgradYdispX + ...
                            locgradYdispY    .* locgradYdispY);

    S11  = lambda_s .* TraE_new .* F11_new + 2 .* miu_s * (F11_new .* E11_new + F12_new * E21_new) ;
    S12  = lambda_s .* TraE_new .* F12_new + 2 .* miu_s * (F11_new .* E12_new + F12_new * E22_new) ;
    S21  = lambda_s .* TraE_new .* F21_new + 2 .* miu_s * (F21_new .* E11_new + F22_new * E21_new) ;
    S22  = lambda_s .* TraE_new .* F22_new + 2 .* miu_s * (F21_new .* E12_new + F22_new * E22_new) ;

    index = 0;
    for i = 1:nen
        for j = 1:nen
            % Galerkin inertia term
            Mij = gW(p)*(N(p,i)*N(p,j));
            Mij = Mij.*solid.dens ;
            Mij = Mij.*volume;
            sA(index+1:index+nElem,p) = Mij;
                    
            % Galerkin convection term

            C1 = gW(p) * (DNDx(:,i) .* DNDx(:,j));
            C2 = gW(p) * (DNDx(:,i) .* DNDy(:,j));
            C3 = gW(p) * (DNDy(:,i) .* DNDx(:,j));
            C4 = gW(p) * (DNDy(:,i) .* DNDy(:,j));

            Aij_1 = ((lambda_s .* (E11 .* dt/2 + E22 .* dt/2 + F11.^2 * dt/2) + ...
                     2 .* nu_s .* (F11.^2 * dt/2 + F12.^2 * dt/4 + E11 .* dt/2)) .* C1 + ...
                     (lambda_s .* (F11 .* F12 * dt/2) + ...
                     2 .* nu_s .* (F11 .* F12 * dt/4 + E21 .* dt/2)) .* C2 + ...
                     (lambda_s .* (F11 .* F21 * dt/2) + ...
                     2 .* nu_s .* (F11 .* F21 * dt/2 + F12 .* F22 .* dt/4)) .* C3 + ...
                     (lambda_s .* (F12 .* F21 * dt/2) + ...
                     2 .* nu_s .* (F11 .* F22 * dt/4)) .* C4) .* volume;
            
            Aij_2 = ((lambda_s .* (F11 .* F21 * dt/2) + ...
                     2 .* nu_s .* (F21 .* F11 * dt/2 + F22 .* F12 * dt/4)) .* C1 + ...
                     (lambda_s .* (F11 .* F22 * dt/2) + ...
                     2 .* nu_s .* (F21 .* F12 * dt/4)) .* C2 + ...
                     (lambda_s .* (E11 .* dt/2 + E22 .* dt/2 + F21.^2 * dt/2) + ...
                     2 .* nu_s .* (F21 .^2 * dt/2 + F22 .^2 * dt/4 + E11 * dt/2)) .* C3 + ...
                     (lambda_s .* (F21 .* F22 * dt/2) + ...
                     2 .* nu_s .* (F21 .* F22 * dt/4 + E21 * dt/2)) .* C4) .* volume;
            
            Aij_3 = ((lambda_s .* (F11 .* F12 * dt/2) + ...
                     2 .* nu_s .* (F11 .* F12 * dt/4 + E12 * dt/2)) .* C1 + ...
                     (lambda_s .* (E11 * dt/2 + E22 * dt/2 + F12.^2 * dt/2) + ...
                     2 .* nu_s .* (F11 .^2 * dt/4 + F12 .^2 * dt/2 + E22 .* dt/2)) .* C2 + ...
                     (lambda_s .* (F11 .* F22 * dt/2) + ...
                     2 .* nu_s .* (F12 .* F21 * dt/4)) .* C3 + ...
                     (lambda_s .* (F12 .* F22 * dt/2) + ...
                     2 .* nu_s .* (F11 .* F21 * dt/4 + F12 .* F22 * dt/2)) .* C4) .* volume;
            
            Aij_4 = ((lambda_s .* (F12 .* F21 * dt/2) + ...
                     2 .* nu_s .* (F22 .* F11 * dt/4)) .* C1 + ...
                     (lambda_s .* (F12 .* F22 * dt/2) + ...
                     2 .* nu_s .* (F21 .* F11 * dt/4 + F22 .* F12 * dt/2)) .* C2 + ...
                     (lambda_s .* (F21 .* F22 * dt/2) + ...
                     2 .* nu_s .* (F21 .* F22 * dt/4 + E12 .* dt/2)) .* C3 + ...
                     (lambda_s .* (E11 * dt/2 + E22 * dt/2 + F22 .^2 * dt/2) + ...
                     2 .* nu_s .* (F21 .^2 * dt/4 + F22 .^2 * dt/2 + E22 * dt/2)) .* C4) .* volume;
            sA1(index+1:index+nElem,p) = Aij_1;
            sA2(index+1:index+nElem,p) = Aij_2;
            sA3(index+1:index+nElem,p) = Aij_3;
            sA4(index+1:index+nElem,p) = Aij_4;

            Aij_5 = gW(p)*N(p,i)*N(p,i)
            Aij_5 = Aij_5 .* solid.dens .* volume;
            sA5(index+1:index+nElem,p) = Aij_5;
            
            Aij_6 = gW(p)*(DNDx(:,i).* S11 + DNDy(:,j) .* S21) *N(p,j);
            Aij_7 = gW(p)*(DNDx(:,i).* S12 + DNDy(:,j) .* S22) *N(p,j);
            sA6(index+1:index+nElem,p) = Aij_6;
            sA7(index+1:index+nElem,p) = Aij_7;

            index = index + nElem;
        end
    end
end
% Summation of all quadrature data
sA = sum(sA,2);
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);
sA7 = sum(sA7,2);
 
% Assemble the matrix
Ms = sparse(iif,jjf,sA,ndof,ndof);  
ZeroF = sparse(ndof,ndof);
        
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);
A6 = sparse(iif,jjf,sA6,ndof,ndof);
A7 = sparse(iif,jjf,sA7,ndof,ndof);

Ms = [Ms ZeroF ;...
      ZeroF Ms ];

Ke = [A1 A2;
      A3 A4]

Src = [A5.*solid.gravity(1) ZeroF
       ZeroF A5.*solid.gravity(2)];

Src1 = [A6]*[conn(1*ndof,1)];
Src2 = [A7]*[conn(1*ndof,1)];
   
MS = (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Ms;

% Left-hand side matrix
LHS = Ke + MS;
        
clear sA Mij sA1 sA2 sA3 Aij_1 Aij_2 Aij_3 sA4 sA5 sA6 sA7 sA8 sA9
clear sA10 sA11 sA12 sA13 sA14 sA15 sA16
clear Aij_4 Aij_5 Aij_6 Aij_7 Aij_8 Aij_9 Aij_10 Aij_11 Aij_12 Aij_13 Aij_14 Aij_15 Aij_16
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16

uAlpha = velSAlpha(:,:,1) ;
uAlpha = [uAlpha(:); zeros(ndof,1)];
uDotAlpha = velSDotAlpha(:,:,1) ;
uDotAlpha = [uDotAlpha(:); zeros(ndof,1)];
if (size(Sol.p,1) == 1 & size(Sol.p,2)>1)
   Sol.p = [Sol.p]';
end

gravVec = [ones(2*ndof,1); zeros(ndof,1)];

% Right-hand side vector
RHS = -(Ms * uDotAlpha(:)) - (Ke)* uAlpha(:) ;
RHS = RHS + (Src)*gravVec ;

                                      
% Solve the linear system
        
% Select the unknown nodal values
freeNodesU = unique([solid.DirichletU; unique(BCCyl(:))]);
freeNodesU = setdiff(1:size(crd,1),[freeNodesU]);
freeNodesV = unique([solid.DirichletV; unique(BCCyl(:))]);
freeNodesV = setdiff(1:size(crd,1),[freeNodesV]);
freeNodesP = setdiff(1:size(crd,1),[]) ;

freeNodes = [freeNodesU';freeNodesV' + size(crd,1); freeNodesP' + 2*size(crd,1)];
        
result = velSAlpha(:,:,1);
result = result(:);
result = [result;Sol.p];
resultDot = velSDotAlpha(:,:,1);
resultDot = resultDot(:);
resultDot = [resultDot; Sol.p];

Increment = LHS(freeNodes,freeNodes)\RHS(freeNodes);
        
% Update the increments
result(freeNodes) = result(freeNodes) + Increment;
resultDot(freeNodes) = resultDot(freeNodes) + (pmc.alphaM/(pmc.gamma*pmc.alpha*solver.dt))*Increment ;

velSAlpha(:,:,1) = reshape(result(1:2*ndof),[],2);
velSDotAlpha(:,:,1) = reshape(resultDot(1:2*ndof),[],2);
Sol.p = result((2*ndof+1):(3*ndof));

% Update the solution
velS = velSPrev + (1/pmc.alpha)*( velSAlpha - velSPrev ) ;
velSDot = velSDotPrev + (1/pmc.alphaM)*( velSDotAlpha - velSDotPrev ) ;
        
NSnormIndicator =  norm(Increment)/norm(result(freeNodes)) ;
fprintf('NS: %e, ', NSnormIndicator);
clear freeNodes1 freeNodes2 freeNodes3
clear result resultDot
        
end

