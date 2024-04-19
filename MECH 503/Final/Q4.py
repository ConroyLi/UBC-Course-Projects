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
    norm_p = np.linalg.norm(p)
    p = p / norm_p
    pb, p2, p3 = p
    
    A11 = C[0,0]*pb**2 + C[0,1]*pb*p2 + C[0,2]*pb*p3 + C[0,1]*p2*pb + C[1,1]*p2**2 + C[1,2]*p2*p3 + C[0,2]*p3*pb + C[1,2]*p3*p2 + C[2,2]*p3**2
    A12 = C[5,0]*pb**2 + C[5,1]*pb*p2 + C[5,1]*pb*p3 + C[5,1]*p2*pb + C[5,1]*p2**2 + C[5,1]*p2*p3 + C[5,1]*p3*pb + C[5,1]*p3*p2 + C[5,1]*p3**2
    A13 = C[0,2]*pb**2 + C[4,1]*pb*p2 + C[4,2]*pb*p3 + C[5,1]*p2*pb + C[2,1]*p2**2 + C[4,3]*p2*p3 + C[4,2]*p3*pb + C[4,3]*p3*p2 + C[4,2]*p3**2
    A22 = C[1,0]*pb**2 + C[1,1]*pb*p2 + C[1,2]*pb*p3 + C[1,1]*p2*pb + C[1,1]*p2**2 + C[1,2]*p2*p3 + C[1,2]*p3*pb + C[1,2]*p3*p2 + C[2,2]*p3**2
    A23 = C[3,0]*pb**2 + C[3,1]*pb*p2 + C[3,2]*pb*p3 + C[3,1]*p2*pb + C[3,1]*p2**2 + C[3,2]*p2*p3 + C[3,2]*p3*pb + C[3,2]*p3*p2 + C[3,2]*p3**2
    A33 = C[2,0]*pb**2 + C[2,1]*pb*p2 + C[2,2]*pb*p3 + C[2,1]*p2*pb + C[2,1]*p2**2 + C[2,2]*p2*p3 + C[2,2]*p3*pb + C[2,2]*p3*p2 + C[2,2]*p3**2
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

s_values = np.linspace(0, 1, 100)
speeds_pb = []
speeds_pp = []

for s in s_values:
    pb = np.array([1-s, s, 0])
    pp = np.array([1-s, 0, s])
    
    A_pb = acoustic_tensor(pb, C)
    A_pp = acoustic_tensor(pp, C)

    speeds_pb.append(wave_speeds(A_pb))
    speeds_pp.append(wave_speeds(A_pp))

speeds_pb = np.array(speeds_pb)
speeds_pp = np.array(speeds_pp)

angles = np.arctan2(s_values, 1 - s_values)  

fig, (ax1, ax2) = plt.subplots(2, 1, subplot_kw={'projection': 'polar'}, figsize=(6, 12))


for i in range(3):  
    ax1.plot(angles, speeds_pb[:, i], label=f"Mode {i+1}")
ax1.set_title("Wave Speeds in Basal Plane")
ax1.set_ylim(0, np.max(speeds_pb) * 1.1)
ax1.legend()


for i in range(3):
    ax2.plot(angles, speeds_pp[:, i], label=f"Mode {i+1}")
ax2.set_title("Wave Speeds in Prismatic Plane")
ax2.set_ylim(0, np.max(speeds_pp) * 1.1)
ax2.legend()

plt.show()