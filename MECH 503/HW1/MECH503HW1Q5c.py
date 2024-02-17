from sympy import symbols, sin, cos, exp, pi, sqrt, diff, Matrix, transpose, evalf

# Define symbols
x1, x2, x3, t = symbols('x1 x2 x3 t')
L = 1
lamda = 1.25

# Identity matrix
I = Matrix.eye(3)

# Define r
r = sqrt(x1**2 + x2**2)

# Define f(x3), phi(x3), y1, y2, y3
f = r + r * sin(pi * x3 / L)
phi = (30 * pi / 180) * exp(2 * x3)
y1 = f * cos(phi) * x1 - f * sin(phi) * x2
y2 = f * sin(phi) * x1 + f * cos(phi) * x2
y3 = lamda * x3

# Displacement
u1 = y1 - x1
u2 = y2 - x2
u3 = y3 - x3

# Gradient of the displacement
du1dx1 = diff(u1, x1)
du1dx2 = diff(u1, x2)
du1dx3 = diff(u1, x3)
du2dx1 = diff(u2, x1)
du2dx2 = diff(u2, x2)
du2dx3 = diff(u2, x3)
du3dx1 = diff(u3, x1)
du3dx2 = diff(u3, x2)
du3dx3 = diff(u3, x3)
gradu = Matrix([[du1dx1, du1dx2, du1dx3],
                [du2dx1, du2dx2, du2dx3],
                [du3dx1, du3dx2, du3dx3]])

# Deformation gradient
myF = I + gradu

# Compute the deformation gradient F
F = Matrix([[diff(y1, x1), diff(y1, x2), diff(y1, x3)],
            [diff(y2, x1), diff(y2, x2), diff(y2, x3)],
            [diff(y3, x1), diff(y3, x2), diff(y3, x3)]])

# Inverse of F
Fm1 = F.inv()

# Green Lagrange strain tensor
EL = 0.5 * (transpose(F) * F - I)

# Euler strain tensor
eE = 0.5 * (I - transpose(Fm1) * Fm1)

# Infinitesimal strain tensor
epsilon = 0.5 * (gradu + transpose(gradu))

# Numerical evaluation of tensors
t0 = 0.05
inf_strain = epsilon.subs({x1: 0.5, x2: 1, x3: 0, t: t0})
Euler_strain = eE.subs({x1: 0.5, x2: 1, x3: 0, t: t0})
Lagrange_strain = EL.subs({x1: 0.5, x2: 1, x3: 0, t: t0})

# Project the tensors over the specific direction m
m = Matrix([1, 1, 0])
epsL = m.transpose() * Lagrange_strain * m
epsE = m.transpose() * Euler_strain * m
epsInf = m.transpose() * inf_strain * m

# Output results
print("Infinitesimal Strain:", epsInf.evalf())
print("Euler Strain:", epsE.evalf())
print("Lagrange Strain:", epsL.evalf())
