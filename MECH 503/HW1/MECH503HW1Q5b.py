
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize

# Initialization
I = np.eye(3)
n0 = np.array([1, 1, 0])
n1 = np.array([0, 1, 1])

xmin, xmax = -1, 1
ymin, ymax = -1, 1
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
    theta = np.linspace(0, 2 * np.pi, 25)
    z = np.linspace(0, height, 25)
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
det_F = np.zeros_like(X)
'''# Deformation of inner surface
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
'''
def calculate_tensors_and_properties(X, Y, Z, lambda_, n0, n1, I):
    y1, y2, y3 = np.zeros_like(X), np.zeros_like(Y), np.zeros_like(Z)
    EL, el, dl, cos_theta = np.zeros_like(X), np.zeros_like(X), np.zeros_like(X), np.zeros_like(X)
    det_F = np.zeros_like(X)
    for i in range(len(X)):
        for j in range(len(X[0])):
            xp = np.array([X[i, j], Y[i, j], Z[i, j]])
            r = np.sqrt(xp[0]**2 + xp[1]**2)
            phi = compute_phi(twist_angle, xp, L)
            f = compute_f(phi, radial_function, R0, xp, L, r0)
           
            F = np.array([[f * np.cos(phi), -f * np.sin(phi), 0],
                          [f * np.sin(phi), f * np.cos(phi), 0],
                          [0, 0, lambda_]])
            F_inv = np.linalg.inv(F)
            yp = F @ xp
            det_F[i, j] = np.linalg.det(F)
            y1[i, j], y2[i, j], y3[i, j] = yp
            EL[i, j] = np.linalg.norm(F.T @ F - I)
            el[i, j] = np.linalg.norm(I - F_inv.T @ F_inv)
            dl[i, j] = np.linalg.norm(F @ n0)
            cos_theta[i, j] = np.dot((F @ n0).T, F @ n1) / (np.linalg.norm(F @ n0) * np.linalg.norm(F @ n1))
            # print(f"Determinant of F at point ({i}, {j}): {det_F}")
    return y1, y2, y3, EL, el, dl, cos_theta

# Compute tensors and properties for inner and outer surfaces
y1, y2, y3, EL, el, dl, cos_theta = calculate_tensors_and_properties(X, Y, Z, lambda_, n0, n1, I)
ym1, ym2, ym3, ELm, elm, DL, cos_THETA = calculate_tensors_and_properties(x, y, z, lambda_, n0, n1, I)

# Normalize the strain tensor values for color mapping
norm = Normalize(vmin=np.min(EL), vmax=np.max(EL))
colors_EL = cm.jet(norm(EL))

norm_m = Normalize(vmin=np.min(ELm), vmax=np.max(ELm))
colors_ELm = cm.jet(norm_m(ELm))

norm_el = Normalize(vmin=np.min(el), vmax=np.max(el))
colors_el = cm.jet(norm_el(el))

norm_elm = Normalize(vmin=np.min(elm), vmax=np.max(elm))
colors_elm = cm.jet(norm_elm(elm))

# Plotting
fig = plt.figure(figsize=(12, 6))
ax1 = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122, projection='3d')

# Plot undeformed and deformed cylinder with strain measures
ax1.plot_surface(X, Y, Z, facecolors=colors_EL, alpha=0.3)
ax1.plot_surface(x, y, z, facecolors=colors_ELm, alpha=0.3)
ax1.set_xlim([xmin, xmax])
ax1.set_ylim([ymin, ymax])
ax1.set_zlim([zmin, zmax])
ax1.set_xlabel('x-axis')
ax1.set_ylabel('y-axis')
ax1.set_zlabel('z-axis')
ax1.set_title('Undeformed Cylinder')

ax2.plot_surface(y1, y2, y3, facecolors=colors_el, alpha=0.3)
ax2.plot_surface(ym1, ym2, ym3, facecolors=colors_elm, alpha=0.3)
ax2.set_xlim([xmin, xmax])
ax2.set_ylim([ymin, ymax])
ax2.set_zlim([zmin, zmax])
ax2.set_xlabel('x-axis')
ax2.set_ylabel('y-axis')
ax2.set_zlabel('z-axis')
ax2.set_title('Deformed Cylinder')


plt.show()
print(det_F[12,12])