N = 1024
dx = 0.01

M = np.arange(1,N+1)
x = -(N-1)*dx/2+(M-1)*dx

one_m3 = np.diag(np.ones(N-3),3)
one_m2 = np.diag(np.ones(N-2),2)
one_m1 = np.diag(np.ones(N-1),1)

one_p3 = np.diag(np.ones(N-3),-3)
one_p2 = np.diag(np.ones(N-2),-2)
one_p1 = np.diag(np.ones(N-1),-1)

up = -2/120*one_m3+18/120*one_m2-3/4*one_m1+3/4*one_p1-18/120*one_p2+2/120*one_p3
u2p = 3/24*one_m3-one_m2+13/8*one_m1-13/8*one_p1+one_p2-3/24*one_p3
u5p = -1/2*one_m3+2*one_m2-5/2*one_m1+5/2*one_p1-2*one_p2+1/2*one_p3


