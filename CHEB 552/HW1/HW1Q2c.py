from sympy import symbols, diff, solve, hessian, Matrix, latex

# Define symbols
F, n, R, P= symbols('F n R P')

# Cost function
F = 14720*(100-P) + (6560- 30.2*P)*(R+1) + 19.5*n*(5000*R - 23*R*P + 5000 - 23*P)**0.5 + 23.2*(5000*R - 23*R*P + 5000 - 23*P)**0.62

# Partial derivatives
DF_R = diff(F,R)
DF_n = diff(F,n)
DF_P = diff(F,P)

# Reported optimum to test
R_asr = 8
n_asr = 55
P_asr = 99

# Test
Result = [DF_R.subs([(R,R_asr),(n,n_asr),(P,P_asr)]),
          DF_n.subs([(R,R_asr),(n,n_asr),(P,P_asr)]),
          DF_P.subs([(R,R_asr),(n,n_asr),(P,P_asr)])]

#is_optimum = True
#if Result == [0,0,0]:
#   print(is_optimum)
#else :
#   is_optimum = False
#print(is_optimum)
print(Result)