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
% ndof = size(crd,1);
% Sol.aleDispPrev = Sol.aleDisp;
Sol.aleDisp = zeros(ndof,2);

% Form the local to global map
iif = zeros(nen^2*nElem,1); 
jjf = zeros(nen^2*nElem,1);
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
eta_Cylx = Sol.dispS(1,1);
eta_Cyly = Sol.dispS(1,2);

BC_n = [BCLeft;BCRight;BCTop;BCBottom;BCCyl]; 

Sol.aleDisp(unique(BCCyl),1) = eta_Cylx; %.*ones(size(unique(BCCyl),1),1);
Sol.aleDisp(unique(BCCyl),2) = eta_Cyly; %.*ones(size(unique(BCCyl),1),1);
Sol.aleDisp(unique(BCTop),1) = eta_xL.*ones(size(unique(BCTop),1),1);
Sol.aleDisp(unique(BCTop),2) = eta_xL.*ones(size(unique(BCTop),1),1);
Sol.aleDisp(unique(BCBottom),1) = eta_x0.*ones(size(unique(BCBottom),1),1);
Sol.aleDisp(unique(BCBottom),2) = eta_x0.*ones(size(unique(BCBottom),1),1);
Sol.aleDisp(unique(BCLeft),1) = eta_y0.*ones(size(unique(BCLeft),1),1);
Sol.aleDisp(unique(BCLeft),2) = eta_y0.*ones(size(unique(BCLeft),1),1);
Sol.aleDisp(unique(BCRight),1) = eta_yL.*ones(size(unique(BCRight),1),1);
Sol.aleDisp(unique(BCRight),2) = eta_yL.*ones(size(unique(BCRight),1),1);
        
% ALE mesh equation
% Localize the data to each element
xxf = zeros(size(cnn));
yyf = zeros(size(cnn));
%dx = zeros(size(cnn));
%dy = zeros(size(cnn));

for i=1:nen
   xxf(:,i) = crd(cnn(:,i),1);
   yyf(:,i) = crd(cnn(:,i),2);
   % dx(:,i) = Sol.dispS(1,1)*ones(nElem,1);
   % dy(:,i) = Sol.dispS(2,1)*ones(nElem,1);
end

% Form element matrix and assemble Galerkin terms
sA1 = zeros(nen^2*nElem,nQuad);
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);
sA4 = zeros(nen^2*nElem,nQuad);
sA5 = zeros(nen^2*nElem,nQuad);

% Quadrature loop
for p = 1:nQuad  
    % Jacobian for each element
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]]; 
    
    % Metric for each element
    volume =abs( J(:,1).*J(:,4) -J(:,2).*J(:,3) );         
    km(:,p) = (max(volume(:)) - min(volume(:)))./volume(:);

    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);
    

    index = 0;
    for i = 1:nen
       for j = 1:nen
           % Galerkin diffusion term
           Aij_1 = gW(p)*(DNDx(:,i).*DNDx(:,j));
           Aij_2 = gW(p)*(DNDy(:,i).*DNDy(:,j));
           Aij_3 = gW(p)*(DNDx(:,i).*DNDy(:,j));
           Aij_4 = gW(p)*(DNDx(:,i)+DNDy(:,i)).*(DNDx(:,j)); 
           Aij_5 = gW(p)*(DNDx(:,i)+DNDy(:,i)).*(DNDy(:,j));
           
           Aij_1 = Aij_1.*volume;
           Aij_2 = Aij_2.*volume;
           Aij_3 = Aij_3.*volume;
           Aij_4 = Aij_4.*volume;
           Aij_5 = Aij_5.*volume;
           
           sA1(index+1:index+nElem,p) = Aij_1.*(1+km(:,p));
           sA2(index+1:index+nElem,p) = Aij_2.*(1+km(:,p));
           sA3(index+1:index+nElem,p) = Aij_3.*(1+km(:,p));
           sA4(index+1:index+nElem,p) = Aij_4.*(1+km(:,p));
           sA5(index+1:index+nElem,p) = Aij_5.*(1+km(:,p));

%            % Galerkin source term
%            Aij_6 = gW(p)*(N(p,i).*N(p,j));
%            Aij_6 = Aij_6.*volume.*f ;
%            sA6(index+1:index+nElem,p) = Aij_6;
           
           index = index + nElem;
       end
    end
end
    
% Summation of all quadrature points
%km = sum(km,2);
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);


% Assemble the matrix      
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A3 = sparse(iif,jjf,sA3,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
A5 = sparse(iif,jjf,sA5,ndof,ndof);
ZeroF = sparse(ndof,ndof);

% Left-hand side matrix
LHS = [2*A1+A2+A4 A3';...
      A3 A1+2*A2+A5  ];

% Right-hand side vector
DVec = LHS*Sol.aleDisp(:);
RHS = [ZeroF; ZeroF] - DVec;

% Select the unknown nodal values
freeNodes = unique([unique(BC_n(:));unique(BCCyl(:))]);
freeNodes1 = setdiff(1:size(crd,1),freeNodes);
% freeNodesY = unique([unique(Bc_n(:));unique(BCCyl(:))]);
% freeNodesY = setdiff(1:size(crd,1),freeNodesY);
freeNodes = [freeNodes1';freeNodes1'+size(crd,1)];
freeNodes1 = [freeNodes1'];

% Solve the linear system
Result = LHS(freeNodes,freeNodes)\RHS(freeNodes);
Result = reshape(Result,[],2);

% Update the ALE displacement and velocity
Sol.aleDisp(freeNodes1,1) = Result(:,1);
Sol.aleDisp(freeNodes1,2) = Result(:,2);
%%
% Sol.aleDisp = Sol.aleDispPrev + (1/pmc.alpha)*( Sol.aleDisp - Sol.aleDispPrev ) ;
%%
% Sol.aleDispy(freeNodes,:) = Result(:,2) ;
Sol.aleVel(freeNodes1,:) = (Sol.aleDisp(freeNodes1,:) - Sol.aleDispPrev(freeNodes1,:))/solver.dt;
% Sol.alevely(freeNodes,:) = Sol.aleDispy(freeNodes,:) - Sol.aleDispprey(freeNodes,:);
end