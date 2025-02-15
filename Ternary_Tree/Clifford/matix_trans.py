import numpy as np

from functools import reduce
n = 8
wires = list(range(n))

P = np.zeros((n,n))
for i in range(n):
    P[wires[i], wires[(i+1)%n]] = 1
    P[wires[(i + 1) % n], wires[(i) % n]] = 1

def T_ij(i, j):
    T = np.eye(n)
    T[i,i] = 0
    T[j,j] = 0
    T[j,i] = 1
    T[i,j] = 1
    return T
T_psn1 = reduce(lambda s,t: s@t, [T_ij(i,(i+1)%n) for i in range(0,n-1,2)], np.eye(n))
T_psn2 = reduce(lambda s,t: s@t, [T_ij(i,(i+1)%n) for i in range(1,n-1,2)], np.eye(n))

print(P)
for i in range(n//2):
    print("i = ", i)
    P = T_psn1.T@P@T_psn1
    P = T_psn2.T@P@T_psn2
    print(P)
