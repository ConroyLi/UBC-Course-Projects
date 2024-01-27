from sympy import symbols, diff, solve, hessian, Matrix, latex

# Define symbols
W, C_p, T_1, P_2, P_1, k, P_3= symbols('W C_p T_1 P_2 P_1 k P_3')

# Work
W = C_p*T_1*((P_2/P_1)**((k-1)/k)+(P_3/P_2)**((k-1)/k)-2)

DW_p2 = diff(W, P_2)
P_2_asr = solve(DW_p2,P_2)
print('W=',latex(W))
print('DW_p2=',latex(DW_p2))
#print(P_2_asr[1])
W_eva = W.subs([(P_2,P_2_asr[1])])
print('P2asr=',latex(P_2_asr))
print('W_eva=',latex(W_eva))