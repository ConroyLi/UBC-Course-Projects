clear all,clc,close all
syms x1 x2 x3 t
L=1;
lamda=1.25;
pi=3.1415926
I = eye(3);
r = sqrt(x1*x1+x2*x2)

f(x3)=r+r*sin(pi*x3/L)
phi(x3)=(30*pi/180)*exp(2*x3)

y1=f(x3)*cos(phi(x3))*x1-f(x3)*sin(phi(x3))*x2;
y2=f(x3)*sin(phi(x3))*x1+f(x3)*cos(phi(x3))*x2;
y3=lamda*x3

%y1 =(r+r*sin(pi*x3/L))*cos((30*pi/180)*exp(2*x3))*x1 - (r+r*sin(pi*x3/L))*sin((30*pi/180)*exp(2*x3))*x2;
%y2 = (r+r*sin(pi*x3/L))*sin((30*pi/180)*exp(2*x3))*x1 + (r+r*sin(pi*x3/L))*cos((30*pi/180)*exp(2*x3))*x2;
%y3 = lamda*x3;
u1 = y1-x1;
u2 = y2-x2;
u3 = y3-x3;
%% Compute the gradient of the displacement
du1dx1 = diff(u1,x1);
du1dx2 = diff(u1,x2);
du1dx3 = diff(u1,x3);
du2dx1 = diff(u2,x1);
du2dx2 = diff(u2,x2);
du2dx3 = diff(u2,x3);
du3dx1 = diff(u3,x1);
du3dx2 = diff(u3,x2);
du3dx3 = diff(u3,x3);
gradu = [du1dx1 du1dx2 du1dx3
du2dx1 du2dx2 du2dx3
du3dx1 du3dx2 du3dx3];
myF = I + gradu;
%% Compute the deformation gradient \dyi/dxj = x \nabla
F11 = diff(y1,x1);
F12 = diff(y1,x2);
F13 = diff(y1,x3);
F21 = diff(y2,x1);
F22 = diff(y2,x2);
F23 = diff(y2,x3);
F31 = diff(y3,x1);
F32 = diff(y3,x2);
F33 = diff(y3,x3);
F = [F11 F12 F13
F21 F22 F23
F31 F32 F33];
Fm1 = inv(F);
%% Green Lagrange strain tensor
EL = 0.5*(transpose(F)*F-I)
%% Euler strain tensor
eE = 0.5*(I-transpose(Fm1)*Fm1)
%% Infinitesimal strain tensor
epsilon = 0.5*(gradu + transpose(gradu))
%% Substitue variables to evaluate numerically the tensors.
t0 = 0.05;
inf_strain = subs(epsilon, [x1,x2,x3,t],[0.5,1,0,t0]);
Euler_strain = subs(eE, [x1,x2,x3,t],[0.5,1,0,t0]);
Lagrange_strain =subs(EL, [x1,x2,x3,t],[0.5,1,0,t0]);
%% Project the tensors over the specific direction m.
m = [1 1 0]'
epsL = m'*Lagrange_strain*m
epsE = m'*Euler_strain*m
epsInf = m'*inf_strain*m