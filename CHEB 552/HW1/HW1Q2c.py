from sympy import symbols, diff, solve, hessian, Matrix, latex

# Define symbols
F, n, R, P= symbols('F n R P')

# Cost function
F = 14720*(100-P) + (6560- 30.2*P)*(R+1) + 19.5*n*(5000*R - 23*R*P + 5000 - 23*P)**0.5 + 23.2*(5000*R - 23*R*P + 5000 - 23*P)**0.62

# Partial derivatives
DF_R = diff(F,R)
DF_n = diff(F,n)
DF_P = diff(F,P)

grad_F = Matrix([[DF_R],[DF_n],[DF_P]])
H_P = Matrix([
    [diff(DF_R, R),diff(DF_R, n),diff(DF_R, P)],
    [diff(DF_n, R),diff(DF_n, n),diff(DF_n, P)],
    [diff(DF_P, R),diff(DF_P, n),diff(DF_P, P)]])
print('F=', latex(F))
print('Grad_f=', latex(grad_F))
print('H_P=', latex(H_P))

# Reported optimum to test
R_asr = 8
n_asr = 55
P_asr = 99

# Test
Result = grad_F.subs({R: R_asr, n: n_asr, P: P_asr})
H_P_eva = H_P.subs({R: R_asr, n: n_asr, P: P_asr})
eigenvalues= H_P_eva.eigenvals()
#is_optimum = True
#if Result == [0,0,0]:
#   print(is_optimum)
#else :
#   is_optimum = False
#print(is_optimum)
print('R=', latex(Result))
print('H_P_eva=',latex(H_P_eva))
print(latex(latex(eigenvalues)))