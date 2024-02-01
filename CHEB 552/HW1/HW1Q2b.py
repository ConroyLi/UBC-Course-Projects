from sympy import symbols, diff, solve, hessian, Matrix, latex

# Define symbols
x_a, x_b, P_p, c_add, S, c_s, Y_p, P= symbols('x_a x_b P_p c_add S c_s Y_p P')

# Cost of addtive
c_add = 20*x_a**2 + 10*x_a + 2

# Cost of steam
c_s = 2*10**(-6)*S**2 + 0.003*S +2

# Yield equation Yp = x_b/x_a
Y_p = 0.0001*S*x_a + 0.001*S + 0.3*x_a + 0.1

# Profit
P = 50*Y_p - c_add - c_s

# Partial derivative
DP_xa = diff(P,x_a)
DP_S = diff(P,S)

SDP_xaxa = diff(DP_xa,x_a)
SDP_xaS = diff(DP_xa,S)
SDP_SS = diff(DP_S,S)

H_P = Matrix([[SDP_xaxa, SDP_xaS],[SDP_xaS, SDP_SS]])

#eigenvalues = H_P.eigenvals
# Checking the definiteness of the Hessian matrices

def check_definiteness(matrix):
    """ Check the definiteness of a matrix """
    eigenvalues = matrix.eigenvals()
    if all(value > 0 for value in eigenvalues):
        return "Positive definite (local minimum)"
    elif all(value < 0 for value in eigenvalues):
        return "Negative definite (local maximum)"
    else:
        return "Indefinite (saddle point)"

# Definiteness of the Hessian matrix for Function a and b
definiteness= check_definiteness(H_P)

print(definiteness)

#print ('P=',latex(P))
print(latex(H_P))