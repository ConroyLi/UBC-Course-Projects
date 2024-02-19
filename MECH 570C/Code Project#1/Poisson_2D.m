%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                               %
%    Solve Poisson equation in two-dimensions   %
%         -\grad^2 u = f,     in \Omega         %
%          u        = 0,     on \Gamma_D        %                            
%        n.\grad u  = 0,     on \Gamma_N        %
%                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       \Omega      : [0,L]x[0,L]               %
%       \Gamma_D    : u = 0 on the four sides   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all

% User-inputs
nElem_x = 10 ;      % Number of elements in X
nElem_y = 10 ;      % Number of elements in Y
nElem = nElem_x * nElem_y ;
f     = 5 ;         % Source value
L     = 1 ;         % Domain length
nen   = 4 ;         % Number of nodes per element
u_x0   = 0 ;        % BC at x=0
u_xL   = 0 ;        % BC at x=L
u_y0   = 0 ;        % BC at y=0
u_yL   = 0 ;        % BC at y=L

% Mesh description: coordinates and connectivity
% Coordinates
[xcrd, ycrd] = ndgrid((0:nElem_x)/(nElem_x),(0:nElem_y)/(nElem_y)); 
crd = [xcrd(:), ycrd(:)];

conn = zeros((nElem_x)*(nElem_y),nen);
for i=1:nElem_y
    for j=1:nElem_x
        ielem = (i-1)*(nElem_x)+j;
        inode = (i-1)*(nElem_x+1)+j;
        % Connectivity
        conn(ielem,:) = [inode inode+1 inode+(nElem_x+2) inode+(nElem_x+1)];
    end
end

% Dirichlet boundary conditions
Bc_y0 = [1:nElem_x+1]';
Bc_xL = [2*(nElem_x+1):nElem_x+1:(nElem_y+1)*(nElem_x+1)]';
Bc_yL = [nElem_y*(nElem_x+1)+nElem_x:-1:nElem_y*(nElem_x+1)+1]';
Bc_x0 = [(nElem_y-1)*(nElem_x+1)+1:-(nElem_x+1):nElem_x+2]';
% BC nodes
Bc_n = [Bc_x0;Bc_xL;Bc_y0;Bc_yL]; 
% BC values
Bc_v = [u_x0.*ones(size(Bc_x0,1),1); u_xL.*ones(size(Bc_xL,1),1) ;
        u_y0.*ones(size(Bc_y0,1),1); u_yL.*ones(size(Bc_yL,1),1) ];     

% Shape functions, gauss points and weights for Numerical integration
% Gauss points
gP = [-1/sqrt(3) -1/sqrt(3)
       1/sqrt(3) -1/sqrt(3)
      -1/sqrt(3)  1/sqrt(3)
       1/sqrt(3) 1/sqrt(3)];       
% Gauss weights
gW = [1 1 1 1];         
% Shape functions
N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;
N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
% Derivative of shape functions
Nx(:,1) = -0.25.*(1-gP(:,2)) ;
Nx(:,2) =  0.25.*(1-gP(:,2)) ;
Nx(:,3) =  0.25.*(1+gP(:,2)) ;
Nx(:,4) = -0.25.*(1+gP(:,2)) ;
Ny(:,1) = -0.25.*(1-gP(:,1)) ;
Ny(:,2) = -0.25.*(1+gP(:,1)) ;
Ny(:,3) =  0.25.*(1+gP(:,1)) ;
Ny(:,4) =  0.25.*(1-gP(:,1)) ;
Nx = Nx' ;
Ny = Ny' ;
% Number of quadrature points
nQuad = size(gW,2);      

% Initialization
ndof = size(crd,1);
Sol.u = zeros(ndof,1);

% Global to local mapping
iif = zeros(nen^2*nElem,1); jjf = zeros(nen^2*nElem,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      iif(index+1:index+nElem) = double(conn(:,i)); 
      jjf(index+1:index+nElem) = double(conn(:,j));  
      index = index + nElem;
   end
end

% Satisfy boundary conditions
Sol.u(Bc_n) = Bc_v ;

%% Poisson equation
xxf = zeros(size(conn));
yyf = zeros(size(conn));
u = zeros(size(conn));

% Localize the variables
for i=1:nen
   xxf(:,i) = crd(conn(:,i),1);
   yyf(:,i) = crd(conn(:,i),2);
   u(:,i) =  Sol.u(conn(:,i),1) ;
end

% Initialize element matrices
sA1 = zeros(nen^2*nElem,nQuad);
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);

% Quadrature loop
for p = 1:nQuad  
    % Jacobian for each element
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]]; 
    % Metric for each element
    volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );         
    
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
    
% Assemble the global matrix
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A3 = sparse(iif,jjf,sA3,ndof,ndof);
    
% Left-hand side matrix and Right-hand side vector
LHS = [A1+A2];
RHS = [A3]*[ones(ndof,1)];
    
% Incorporate Dirichlet condition on RHS
DVec = LHS*Sol.u ;
RHS = RHS - DVec ;

% Get the unknown nodal information
freeNodes = unique([Bc_n]);
freeNodes = setdiff(1:size(crd,1),[freeNodes]);
    
% Solve the linear system
Result = LHS(freeNodes,freeNodes)\RHS(freeNodes);

% Update the solution
Sol.u(freeNodes,:) = Result ;

% Post-processing
figure(1)
for i=1:nElem_y+1
    posi = [(i-1)*(nElem_x+1)+1:i*(nElem_x+1)];
    xx(i,:) = crd(posi,1)';
    yy(i,:) = crd(posi,2)';
    Z(i,:) = Sol.u(posi)';
end
surf(xx,yy,Z);
