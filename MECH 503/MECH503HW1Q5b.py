'''
from sympy import zeros, Matrix, latex, sin, cos, pi, exp, sqrt, evalf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
# Define symbols
# x1, x2, x3, y1, y2, y3= symbols('x1 x2 x3 y1 y2 y3 ')

landa = 1.25
r_0 = 0.25
r_1 = 0.75
dr = r_1-r_0
L = 1
N = 2
x1 = [0]
x2 = [0]
x3 = [0]
y1 = [0]
y2 = [0]
y3 = [0]
print(y1)
# Reference config
N = 5
x1[0] = r_0
x2[0] = r_1
x3[0] = 0
for i in range(1,N+1):
    # x1[i+1] = x1[i] + r_0/N
    x1.append(x1[i-1] + dr/N)
    for j in range(1,N+1):
        x2.append(x2[i-1] + dr/N)
        # x2[i+1] = x2[i] + r_0/N
        for k in range(1,N+1):
            x3.append(x3[i-1] + L/N)
            # x3[i+1] = x3[i] + L/N

            r = sqrt(x1[i]**2 + x2[j]**2)
            f = r + r*sin(pi*x3[k]/L)
            phi = (pi/6)*(exp(2*x3[k]))
            
            # Deformation gradient
            F = Matrix([
            [f*cos(phi), -f*sin(phi), 0],
            [f*sin(phi),  f*cos(phi), 0],
            [0,           0     , landa]])

            y = F*Matrix([x1[i],x2[j],x3[k]])
            # print(y)
            dy1 = y[0,0].evalf()
            dy2 = y[1,0].evalf()
            dy3 = y[2,0].evalf()
            y1.append(dy1)
            y2.append(dy1)
            y3.append(dy1)
            #print(x1,x2,x3,y1,y2,y3)
            
            print(y1)
'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Initialization
I = np.eye(3)
n0 = np.array([1, 1, 0])
n1 = np.array([0, 1, 1])

xmin, xmax = -2, 2
ymin, ymax = -2, 2
zmin, zmax = 0, 1.5

radial_function = 'sinusoidal'  # Options: 'linear', 'quadratic', 'sinusoidal'
twist_angle = 'exponential'  # Options: 'constant', 'sinusoidal', 'quadratic', 'exponential'

a, b = 0.25, 0.75  # Radii
r0 = a
R0 = b
L = 1  # Height
lambda_ = 1.25  # Compression/stretching value

# Function to create cylinder data
def create_cylinder(radius, height):
    theta = np.linspace(0, 2 * np.pi, 50)
    z = np.linspace(0, height, 50)
    theta, z = np.meshgrid(theta, z)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    return x, y, z

# Create cylinders
X, Y, Z = create_cylinder(a, L)
x, y, z = create_cylinder(b, L)

# Helper functions for f and phi
def compute_f(phi, radial_function, R0, xp, L, r0):
    if radial_function == 'linear':
        return R0 + 2 * xp[2]
    elif radial_function == 'quadratic':
        return R0 + xp[2]**2
    elif radial_function == 'sinusoidal':
        return r0 + r0 * np.sin(np.pi * xp[2] / L)
    else:
        raise ValueError("Invalid radial function")

def compute_phi(twist_angle, xp, L):
    if twist_angle == 'constant':
        return 30 * np.pi / 180
    elif twist_angle == 'sinusoidal':
        return 30 * np.sin(np.pi * xp[2] / L) * np.pi / 180
    elif twist_angle == 'quadratic':
        return 30 * xp[2]**2 * np.pi / 180
    elif twist_angle == 'exponential':
        return 30 * np.exp(2 * xp[2]) * np.pi / 180
    else:
        raise ValueError("Invalid twist angle")

# Initialize arrays for deformed coordinates
y1, y2, y3 = np.zeros_like(X), np.zeros_like(Y), np.zeros_like(Z)
ym1, ym2, ym3 = np.zeros_like(x), np.zeros_like(y), np.zeros_like(z)

# Deformation of inner surface
for i in range(len(X)):
    for j in range(len(X[0])):
        xp = np.array([X[i, j], Y[i, j], Z[i, j]])
        phi = compute_phi(twist_angle, xp, L)
        f = compute_f(phi,radial_function, R0, xp, L, r0)
       
        F = np.array([[f * np.cos(phi), -f * np.sin(phi), 0],
                      [f * np.sin(phi), f * np.cos(phi), 0],
                      [0, 0, lambda_]])
        yp = F @ xp
        y1[i, j], y2[i, j], y3[i, j] = yp

# Deformation of outer surface
for i in range(len(x)):
    for j in range(len(x[0])):
        xp = np.array([x[i, j], y[i, j], z[i, j]])
        phi = compute_phi(twist_angle, xp, L)
        f = compute_f(phi, radial_function, R0, xp, L, r0)
      
        F = np.array([[f * np.cos(phi), -f * np.sin(phi), 0],
                      [f * np.sin(phi), f * np.cos(phi), 0],
                      [0, 0, lambda_]])
        yp = F @ xp
        ym1[i, j], ym2[i, j], ym3[i, j] = yp

# Plotting
fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

# Plot undeformed cylinder
ax1.plot_surface(X, Y, Z, color='blue', alpha=0.5)
ax1.plot_surface(x, y, z, color='red', alpha=0.5)
ax1.set_xlim([xmin, xmax])
ax1.set_ylim([ymin, ymax])
ax1.set_zlim([zmin, zmax])
ax1.set_xlabel('x-axis')
ax1.set_ylabel('y-axis')
ax1.set_zlabel('z-axis')
ax1.set_title('Undeformed Cylinder')

# Plot deformed cylinder
ax2.plot_surface(y1, y2, y3, color='blue', alpha=0.5)
ax2.plot_surface(ym1, ym2, ym3, color='red', alpha=0.5)
ax2.set_xlim([xmin, xmax])
ax2.set_ylim([ymin, ymax])
ax2.set_zlim([zmin, zmax])
ax2.set_xlabel('x-axis')
ax2.set_ylabel('y-axis')
ax2.set_zlabel('z-axis')
ax2.set_title('Deformed Cylinder')

plt.show()
