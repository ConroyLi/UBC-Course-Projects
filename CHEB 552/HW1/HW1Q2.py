from sympy import symbols, diff, solve, hessian, Matrix, latex

# Define symbols
alpha, c, c_c, c_i, c_x, i, n, p ,q, T= symbols('alpha c c_c c_i c_x i n p q T')

# Linear sum of cost
temp_var= 52.47*q*360 # Temporery value of the denominator
c = c_c + c_i + c_x + (2.09*10**4/360)*T**(-0.3017) + (1.064*10**6*(alpha*T**0.4925/q))/temp_var
+ (4.242*10**4*alpha*T**0.7952+1.813*i**p*(n*T+1.2*q)**0.861)/temp_var
+ (4.25*10**3*alpha*(n*T+1.2*q))/temp_var + ((5.042*10**3*q**(-0.1899))+0.1049*q**0.671)/360

print(latex(c))