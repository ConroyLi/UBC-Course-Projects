"""
import math
import numpy as np

n=5
m= 8
# make a random (m,n) matrix
A= np.reshape( np.random.randint(0, 9, size= m*n), (m, n)).copy()
A.resize((2,2))
print("A=",A)
#SVD 
U,S,Vh = np.linalg.svd(A) 
S_m=np.diag(S)
S_m.resize((U.shape[1], Vh.shape[0]))
A_svd= U@S_m@Vh
#print(A_svd)
#print(np.allclose(A,A_svd))
print(S)
print(S_m)
"""
import numpy as np

# Specify the size of the matrix
n = 5
m = 8

# Create a random (m, n) matrix
matrix_A = np.random.randint(0, 9, size=(m, n))

# Perform SVD
U, singular_values, Vh = np.linalg.svd(matrix_A)

# Construct the singular value matrix
S_m = np.zeros((U.shape[1], Vh.shape[0]))
np.fill_diagonal(S_m, singular_values)

# Reconstruct matrix from SVD
A_svd = U @ S_m @ Vh

# Output
print("Original Matrix A:\n", matrix_A)
print("Reconstructed Matrix from SVD:\n", A_svd)
print("Singular Values:\n", singular_values)
print("Diagonal Matrix of Singular Values:\n", S_m)
print("SVD Reconstruction Close to Original:", np.allclose(matrix_A, A_svd))
