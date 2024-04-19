import numpy as np
import matplotlib.pyplot as plt

C11 = C22 = 63.5  # GPa
C33 = 66.5
C12 = 25.9
C13 = C23 = 21.7
C44 = 18.4
C66 = (C11 - C12) / 2  

C = np.array([
    [C11, C12, C13, 0,   0,   0],
    [C12, C11, C13, 0,   0,   0],
    [C13, C13, C33, 0,   0,   0],
    [0,   0,   0,   C44, 0,   0],
    [0,   0,   0,   0,   C44, 0],
    [0,   0,   0,   0,   0,   C66]
])

def acoustic_tensor(p, C):
    # Unpack direction components
    norm_p = np.linalg.norm(p)
    p = p / norm_p if norm_p != 0 else p
    p1, p2, p3 = p
    
    # Calculate acoustic tensor components using provided relations
    A11 = C[0,0]*p1**2 + C[0,1]*p1*p2 + C[0,2]*p1*p3 + C[0,1]*p2*p1 + C[1,1]*p2**2 + C[1,2]*p2*p3 + C[0,2]*p3*p1 + C[1,2]*p3*p2 + C[2,2]*p3**2
    A12 = C[5,0]*p1**2 + C[5,1]*p1*p2 + C[5,1]*p1*p3 + C[5,1]*p2*p1 + C[5,1]*p2**2 + C[5,1]*p2*p3 + C[5,1]*p3*p1 + C[5,1]*p3*p2 + C[5,1]*p3**2
    A13 = C[0,2]*p1**2 + C[4,1]*p1*p2 + C[4,2]*p1*p3 + C[5,1]*p2*p1 + C[2,1]*p2**2 + C[4,3]*p2*p3 + C[4,2]*p3*p1 + C[4,3]*p3*p2 + C[4,2]*p3**2
    A22 = C[1,0]*p1**2 + C[1,1]*p1*p2 + C[1,2]*p1*p3 + C[1,1]*p2*p1 + C[1,1]*p2**2 + C[1,2]*p2*p3 + C[1,2]*p3*p1 + C[1,2]*p3*p2 + C[2,2]*p3**2
    A23 = C[3,0]*p1**2 + C[3,1]*p1*p2 + C[3,2]*p1*p3 + C[3,1]*p2*p1 + C[3,1]*p2**2 + C[3,2]*p2*p3 + C[3,2]*p3*p1 + C[3,2]*p3*p2 + C[3,2]*p3**2
    A33 = C[2,0]*p1**2 + C[2,1]*p1*p2 + C[2,2]*p1*p3 + C[2,1]*p2*p1 + C[2,1]*p2**2 + C[2,2]*p2*p3 + C[2,2]*p3*p1 + C[2,2]*p3*p2 + C[2,2]*p3**2
    A = np.array([
        [A11, A12, A13],
        [A12, A22, A23],
        [A13, A23, A33]
    ])
    return A
def wave_speeds(A):
    
    eigenvalues = np.linalg.eigvalsh(A)
    print(eigenvalues)
    return np.sqrt(np.maximum(eigenvalues, 0))

# Values of s from 0 to 1
s_values = np.linspace(0, 1, 100)
speeds_p5 = []
speeds_p6 = []

for s in s_values:
    # Direction in the basal plane
    p5 = np.array([1-s, s, 0])
    # Direction in the prismatic plane
    p6 = np.array([1-s, 0, s])
    
    # Compute the acoustic tensors
    A_p5 = acoustic_tensor(p5, C)
    A_p6 = acoustic_tensor(p6, C)
    
    # Compute wave speeds
    speeds_p5.append(wave_speeds(A_p5))
    speeds_p6.append(wave_speeds(A_p6))

speeds_p5 = np.array(speeds_p5)
speeds_p6 = np.array(speeds_p6)

# Prepare polar plots
angles = np.arctan2(s_values, 1 - s_values)  # Angle for polar plot

fig, (ax1, ax2) = plt.subplots(2, 1, subplot_kw={'projection': 'polar'}, figsize=(6, 12))

# Plotting for basal plane
for i in range(3):  # There can be up to 3 wave speeds (longitudinal and two shear)
    ax1.plot(angles, speeds_p5[:, i], label=f"Mode {i+1}")
ax1.set_title("Wave Speeds in Basal Plane")
ax1.set_ylim(0, np.max(speeds_p5) * 1.1)
ax1.legend()

# Plotting for prismatic plane
for i in range(3):
    ax2.plot(angles, speeds_p6[:, i], label=f"Mode {i+1}")
ax2.set_title("Wave Speeds in Prismatic Plane")
ax2.set_ylim(0, np.max(speeds_p6) * 1.1)
ax2.legend()

plt.show()