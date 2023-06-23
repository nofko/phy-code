

def zf(E):

    M = B(u[n][:-1],ux,ut,E[0]+E[1]*1j)
    tr = np.trace(M)
    y = (tr.real)**2+(tr.imag)**2
    
    return [y] * len(E)

solution = optimize.root(zf, [-0.12, 0.21])

print(solution)
