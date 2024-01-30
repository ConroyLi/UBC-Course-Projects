from sympy import symbols, diff, solve, hessian, Matrix, latex

# Define symbols
u, r, e, sigma = symbols('u r e sigma')

# Potential Function
u = 4*e*((sigma/r)**12 - (sigma/r)**6)

# First derivative of u with respect to r
du_r = diff(u, r)

# Solve du_r = 0
r_asr = solve(du_r, r)

# Second derivative of u with respect to r
sdu_r = diff(du_r, r)

# Number of solutions
l = len(r_asr)

# print("Solutions:", r_asr, "Number of solutions:", l)
print('u=',latex(u))
print('r_asr=',latex(r_asr[0:2]))
print('du_r=', latex(du_r))
print('sdu_r=', latex(sdu_r))
# Loop over each solution
'''for i in range(l): # Test first two solutions
    
    sdu_r_eva = sdu_r.subs([(r, r_asr[i].subs(sigma,1)) , (sigma, 1)])
    # Check if the second derivative is a real number
    if sdu_r_eva.is_real:
        if sdu_r_eva > 0:
            print(sdu_r_eva, 'is a local minimum at r =', r_asr[i])
        elif sdu_r_eva < 0:
            print(sdu_r_eva, 'is a local maximum at r =', r_asr[i])
        else:
            print(sdu_r_eva, 'is a saddle point at r =', r_asr[i])
    else:
        print(sdu_r_eva,'Complex second derivative at r =', r_asr[i])
'''
sdu_r_eva_1 = sdu_r.subs(r, r_asr[0])
sdu_r_eva_2 = sdu_r.subs(r, r_asr[1])

u_asr = u.subs(r, r_asr[0])

# print('1=',latex(sdu_r_eva_1))
# print('2=',latex(sdu_r_eva_2))
print('2=',latex(u_asr))