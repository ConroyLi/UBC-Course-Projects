from sympy import symbols, diff, solve, hessian, Matrix, latex, eye, sqrt, Trace

# Define symbols
x1, x2, x3, u, x0, sr, e= symbols('x1 x2 x3 u x0 sr e')

# u
u = Matrix([[x1**3+x1*x2**2+x3**2],
            [x3**2-x2**3],
            [2*x1*x2*x3-x1**3]])
# print('u=',latex(u))

# Grad u and its transpose
grad_u = Matrix([[diff(u[0],x1),diff(u[0],x2),diff(u[0],x3)],
          [diff(u[1],x1),diff(u[1],x2),diff(u[1],x3)],
          [diff(u[2],x1),diff(u[2],x2),diff(u[2],x3)]]) 
# print('grad_u=',latex(grad_u))

grad_uT = grad_u.T
# print('grad_uT=',latex(grad_uT))

# Infinitisimal Strain Tensor
e = 0.5*(grad_u + grad_uT)
e_eva = e.subs({x1:1,x2:0,x3:1})
# print('e_eva=',latex(e_eva))

# Vorticity Tensor
w = 0.5*(grad_u - grad_uT)
w_eva = w.subs({x1:1,x2:0,x3:1})
# print('w=',latex(w))
# print('w_eva=',latex(w_eva))

# Length Change
F = grad_u + eye(3)
# print(latex(F))
x0 = Matrix([[1/sqrt(2)],[0],[1/sqrt(2)]])
sr = F.subs({x1:1,x2:0,x3:1})*x0
# print(latex(sr))
SR = sqrt(sr[0]**2+sr[1]**2+sr[2]**2)
# print(latex(SR-1))

# Volume Change
dV = Trace(e_eva).simplify()
# print('dV=',latex(dV))

# Deviatoric infinitisimal strain tensor is determined to be: