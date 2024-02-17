from sympy import symbols, diff, solve, hessian, Matrix, latex, eye, sqrt, Trace

# Define symbols
EL, F, f= symbols('El F f')

EL = Matrix([[2, 0.1, 0],
             [0.1, 1, 0],
             [0, 0, 0]])
print(latex(EL))
f= [1,2,0]
l0 = sqrt(f[0]**2+f[1]**2+f[2]**2)
n = Matrix([[1/l0],[2/l0],[0]])

eL = n.T*EL*n
# print(eL)
l = l0*eL[0]
print(l)