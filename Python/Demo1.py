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