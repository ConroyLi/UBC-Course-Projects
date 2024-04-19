from sympy import diff, symbols, latex, simplify,integrate,lambdify
import matplotlib.pyplot as plt
import numpy as np
phi, h, L, x, y, q= symbols('phi h l x y q')

phi = q/(20*h**3) * (-20*y**3*(L**2-x**2) - 4*y**5 -15*h**2*x**2*y + 2*h**2*y**3 - 5*h**3*x**2)

sigma_x = diff(diff(phi,y),y)
sigma_y = diff(diff(phi,x),x)
sigma_xy = -diff(diff(phi,x),y)

BC11 = sigma_xy.subs(y,-h/2)
BC12 = sigma_xy.subs(y,h/2)
BC2 = sigma_y.subs(y,-h/2)
BC3 = sigma_y.subs(y,h/2)

BC4 = integrate(sigma_x,(y,-h/2,h/2))
BC41 = BC4.subs(x,-L)
BC42 = BC4.subs(x,L)

BC5 = integrate(sigma_x*y,(y,-h/2,h/2))
BC51 = BC5.subs(x,-L)
BC52 = BC5.subs(x,L)


BC6 = integrate(sigma_xy,(y,-h/2,h/2))
BC61 = BC6.subs(x,-L)
BC62 = BC6.subs(x,L)

print('BC11=',latex(simplify(BC11)))
print('BC12=',latex(simplify(BC12)))
print('BC2=',latex(simplify(BC2)))
print('BC3=',latex(simplify(BC3)))
print('BC41=',latex(simplify(BC41)))
print('BC42=',latex(simplify(BC42)))
print('BC51=',latex(simplify(BC51)))
print('BC52=',latex(simplify(BC52)))
print('BC61=',latex(simplify(BC61)))
print('BC62=',latex(simplify(BC62)))

sigma_xx = sigma_x.subs(x,L)
print('sigma_xx=',latex(simplify(sigma_xx)))

L_value = 2 
h_value = 1   
q_value = 1 
x_values = np.linspace(-L_value, L_value, 400)
bending = integrate(sigma_x*y,(y,-h/2,h/2))
#sigma_x_numerical = lambdify((x, y, h, q, L), sigma_x, 'numpy')
sigma_xy_numerical = lambdify((x, y, h, q, L), sigma_xy, 'numpy')
bending_numerical = lambdify((x, h, q, L), bending, 'numpy')

y_mid = 0


#sigma_x_values_mid_y = sigma_x_numerical(x_values, y_mid, h_value, q_value, L_value)
sigma_xy_values_mid_y = sigma_xy_numerical(x_values, y_mid, h_value, q_value, L_value)
bending_values = bending_numerical(x_values, h_value, q_value, L_value)

plt.figure(figsize=(10, 5))
#plt.plot(x_values, sigma_x_values_mid_y, label=r'$\sigma_x$ at $y = h/2$')
plt.plot(x_values, sigma_xy_values_mid_y, label=r'Shear Stress at $y = 0$')
plt.plot(x_values, bending_values, label=r'Moment at $y = 0$')
plt.title('Shear Stress and Moment Distribution along x-direction at $y=0$')
plt.xlabel('x (Position along the beam)')
plt.ylabel('Normal Stress $\sigma_x$')
plt.axhline(0, color='gray', lw=0.5)  # Add a line at sigma_x = 0 for reference
plt.legend()
plt.grid(True)
plt.show()
