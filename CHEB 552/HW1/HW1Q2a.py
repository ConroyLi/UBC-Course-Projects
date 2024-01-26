from sympy import symbols, diff, solve, hessian, Matrix, latex

# Define symbols
alpha, c, c_c, c_i, c_x, i, n, p ,q, T= symbols('alpha c c_c c_i c_x i n p q T')

# Values of each variable
alpha = 0.2
c_c = 12.5
c_i = 0.5
c_x = 0.9
i = 0.1
n = 2
p = 7000

# Linear sum of cost
# temp_var= # Temporery value of the denominator
# q = q/6.29
c = c_c + c_i + c_x + (2.09*10**4/360)*T**(-0.3017) + (1.064*10**6*(alpha*T**0.4925))/(52.47*q*6.29*360) + (4.242*10**4*alpha*T**0.7952+1.813*i**p*(n*T+1.2*q*6.29)**0.861)/(52.47*q*6.29*360) + (4.25*10**3*alpha*(n*T+1.2*q*6.29))/(52.47*q*6.29*360) + ((5.042*10**3*(q*6.29)**(-0.1899))+0.1049*(q*6.29)**0.671)/360

# print('c=', latex(c))

# Partial derivative of cost function w.r.t q and T
PD_c_q = diff(c,q)
PD_c_T = diff(c,T)
g_c = Matrix([PD_c_q, PD_c_T])
#print('PD_c_q=',latex(PD_c_q),'PD_c_T=', latex(PD_c_T))
# print('PD_c_q=',PD_c_q,'PD_c_T=', PD_c_T)

# Stationary point of cos function
#SP_c = solve((PD_c_q,PD_c_T),(q,T))
#print('stationary points are', latex(SP_c))

# Second partial derivative of cost function
SPD_c_qq = diff(PD_c_q,q)
SPD_c_qT = diff(PD_c_q,T)
SPD_c_TT = diff(PD_c_T,T)
 
# Hessian matrix of cost function
H_c = Matrix([[SPD_c_qq, SPD_c_qT],
              [SPD_c_qT, SPD_c_TT]])
#print('H_c=', latex(H_c))

# Evaluate Hessian matrix 
j = 1
x = Matrix([100, 100]) # initial guess of q and T
while (j>0):
    
    H_c_eva = H_c.subs({q:x[0,0],T:x[1,0]})
    g_c_eva = g_c.subs({q:x[0,0],T:x[1,0]})
    dx = H_c_eva.inv('CH') * g_c_eva
    # 
    x = x - dx
    
    if -dx[0,0]/x[0,0]< 10**(-6):
        print(x)
        break
    j += 1

print('q =',latex(x[0,0]*6.29), 'T=', latex(x[1,0]))

c_min = c.subs([(q,x[0,0]),(T,x[1,0])])

print(c_min)
