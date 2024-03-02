%% Q2
F = [1 0 2; 0 1 -2; -2 2 1];
C = F'*F
m = [1 ; 0; 0];
n = [0 ; 1; 0];
costheta = (m' * C * n)/(sqrt(m' * C * m)*sqrt(n' * C * n))

n_p = [1\sqrt(2);1\sqrt(2);0];
l_p = sqrt(2);
dl = l_p*sqrt(n_p' * C * n_p)

EL = 0.5*(C-eye(3))
EE = 0.5*(eye(3)-inv(C))

[U, S, V] = svd(F);
U_stretch = V*S*V'
R = F*inv(U_stretch)
RT = R'
R_inv = inv(R)
%det = det(R)

[ps,pd] = eig(U_stretch)

%% Q2
S = [-50 50 70; 50 100 0;70 0 70]; %\sigma
c = [1 0 0.5];
n = c'/sqrt(1+0.5^2);
t = S*n
eig(S)

%% Q5
l = [0 1 0];
m = [0 0 1]';
lm = l'*m'
gamma = 1;
F5 = eye(3)+gamma*lm

m_d = F5*m