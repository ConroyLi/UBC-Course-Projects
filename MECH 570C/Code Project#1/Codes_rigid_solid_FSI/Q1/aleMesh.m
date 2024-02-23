%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%                      2D ALE Mesh Equation Solver                        %
%                                                                         % 
%              \nabla . \sigma^m = 0,    in \Omega^s,                     %                      
%              \sigma = grad \eta + grad \eta^T + (div \eta)I,            %
%                                                                         %
%    where \sigma is the stress experienced by the ALE mesh due to        %
%    strain induced by structural movement, \eta is the mesh disp-        %
%    lacement for the fluid nodes and I is the identity tensor.           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sol] = aleMesh(Sol, solver, solid, BCCyl, BCTop, BCBottom,...
                          BCLeft, BCRight, pmc, cnn, crd, elemType, nen, ...
                          ndof, nElem)

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
 
% Initialize the displacement
ndof = size(crd,1);
Sol.u = zeros(ndof,1);
% Form the local to global map
iif = zeros(nen^2*nElem,1); jjf = zeros(nen^2*nElem,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      iif(index+1:index+nElem) = double(cnn(:,i)); % cnn
      jjf(index+1:index+nElem) = double(cnn(:,j));  
      index = index + nElem;
   end
end

% Satisfy boundary conditions
eta_x0 = 0;
eta_xL = 0;
eta_y0 = 0;
eta_yL = 0;
%eta_Cylx = Sol.dispS(unique(BCCyl),1);
%eta_Cyly = Sol.dispS(unique(BCCyl),2);

Sol.aleDisp(unique(BCCyl),1) = eta_Cylx;
Sol.aleDisp(unique(BCCyl),2) = eta_Cyly;
Sol.aleDisp(unique(BCTop),1) = eta_xL.*ones(size(BCTop),1);
Sol.aleDisp(unique(BCTop),2) = eta_xL.*ones(size(BCTop),2);
Sol.aleDisp(unique(BCBottom),1) = eta_x0.*ones(size(BCBottom),1);
Sol.aleDisp(unique(BCBottom),2) = eta_x0.*ones(size(BCBottom),2);
Sol.aleDisp(unique(BCLeft),1) = eta_y0.*ones(size(BCLeft),1);
Sol.aleDisp(unique(BCLeft),2) = eta_y0.*ones(size(BCLeft),2);
Sol.aleDisp(unique(BCRight),1) = eta_yL.*ones(size(BCRight),1);
Sol.aleDisp(unique(BCRight),2) = eta_yL.*ones(size(BCRight),2);
        
% ALE mesh equation
% Localize the data to each element
xxf = zeros(size(cnn));
yyf = zeros(size(cnn));
u = zeros(size(cnn));
for i=1:nen
   xxf(:,i) = crd(cnn(:,i),1);
   yyf(:,i) = crd(cnn(:,i),2);
   u(:,i) =  Sol.aleDisp(cnn(:,i),1) ;
end

% Form element matrix and assemble Galerkin terms
sA1 = zeros(nen^2*nElem,nQuad);
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);

% Quadrature loop
for p = 1:nQuad  
    % Jacobian for each element
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]]; 
    % Metric for each element
    volume =abs( J(:,1).*J(:,4) -J(:,2).*J(:,3) );         
    
    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    
    index = 0;
    for i = 1:nen
       for j = 1:nen
           % Galerkin diffusion term
           Aij_1 = gW(p)*(DNDx(:,i).*DNDx(:,j));
           Aij_2 = gW(p)*(DNDy(:,i).*DNDy(:,j));
           Aij_1 = Aij_1.*volume;
           Aij_2 = Aij_2.*volume;
           sA1(index+1:index+nElem,p) = Aij_1;
           sA2(index+1:index+nElem,p) = Aij_2;
           
           % Galerkin source term
           Aij_3 = gW(p)*(N(p,i).*N(p,j));
           Aij_3 = Aij_3.*volume.*f ;
           sA3(index+1:index+nElem,p) = Aij_3;
           
           index = index + nElem;
       end
    end
end
    
% Summation of all quadrature points
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);

% Summation of all quadrature data
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
 
% Assemble the matrix      
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A3 = sparse(iif,jjf,sA3,ndof,ndof);

% Left-hand side matrix
LHS = [A1+A2];

% Right-hand side vector
RHS = [A3]*[ones(ndof,1)];

% Select the unknown nodal values
freeNodes = unique([Bc_n]);
freeNodes = setdiff(1:size(crd,1),[freeNodes]);

% Solve the linear system
Result = LHS(freeNodes,freeNodes)\RHS(freeNodes);

% Update the ALE displacement and velocity
Sol.u(freeNodes,:) = Result ;
Sol.aleVel = Sol.aleDisp/solver.dt;
end
