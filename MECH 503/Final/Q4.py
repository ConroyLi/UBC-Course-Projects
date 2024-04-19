import numpy as np
import matplotlib.pyplot as plt

C11 = C22 = 63.5  # GPa
C33 = 66.5
C12 = 25.9
C13 = C23 = 21.7
C44 = C55 = 18.4
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
  
    p_ext = np.array([p[0], p[1], p[2], p[1]*p[2], p[0]*p[2], p[0]*p[1]])
    
    Q = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            Q[i, j] = np.sum(C[i*2:i*2+2, j*2:j*2+2] * np.outer(p_ext[i*2:i*2+2], p_ext[j*2:j*2+2]))
    return Q

def wave_speeds(Q):
    
    eigenvalues = np.linalg.eigvalsh(Q)
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
    Q_p5 = acoustic_tensor(p5, C)
    Q_p6 = acoustic_tensor(p6, C)
    
    # Compute wave speeds
    speeds_p5.append(wave_speeds(Q_p5))
    speeds_p6.append(wave_speeds(Q_p6))

speeds_p5 = np.array(speeds_p5)
speeds_p6 = np.array(speeds_p6)

# Prepare polar plots
angles = np.arctan2(s_values, 1 - s_values)  # Angle for polar plot

fig, (ax1, ax2) = plt.subplots(1, 2, subplot_kw={'projection': 'polar'}, figsize=(12, 6))

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