## Finding the  Estar  -  possibly can be improved, the issue is the case when we are missing a mode

Estar = (E[np.where(np.diff(S)==1)[0]]+E[np.where(np.diff(S)==1)[0]+1])/2
Estar_idx = np.where(np.diff(S)==1)[0]


plt.figure()
plt.title("$M_{11}$ element of the matrix")
plt.plot(E,M11)
plt.plot(E,S)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("$M_{11}$")

plt.scatter(E[Estar_idx],np.zeros_like(Estar))



## We compute the M12 for each of the Estar

M12star = np.zeros(len(Estar))


def M12_finder(E):
    Mt, sst = B(lam*u+E)
    
    return Mt[0,1]

def trace_finder_plus(E):
    Mt, sst = B(lam*u+E)

    tr = np.trace(Mt)/2
    
    if np.abs(tr)<=1:
        return tr-1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))-1
    
def trace_finder_minus(E):
    Mt, sst = B(lam*u+E)

    tr = np.trace(Mt)/2
    
    if np.abs(tr)<=1:
        return tr+1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))+1

for i in range(len(Estar)):
    
    M12star[i] = M12_finder(Estar[i])


mus = [E_min]

for i in range(len(Estar)-1):

    mus.append(bisect(M12_finder,Estar[i],Estar[i+1]))

mus=np.array(mus)

plt.figure()
plt.title("$M_{12}$ element of the matrix")
plt.plot(E,M12)
plt.plot(E,S)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axhline(0,0,1,ls="--",color="tab:gray")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("$M_{12}$")
plt.scatter(mus,np.zeros_like(mus))


Eigs_plus = []
Eigs_minus = []


for i in range(len(mus)-1):

    Eigs_plus.append(bisect(trace_finder_plus,mus[i],mus[i+1],xtol=1e-12))
    Eigs_minus.append(bisect(trace_finder_minus,mus[i],mus[i+1],xtol=1e-12))


plt.figure()
plt.title("Trace")
plt.plot(E,trace)
plt.plot(E,S)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.xlabel("E")
plt.ylabel("Trace")

plt.scatter(Eigs_plus,np.ones_like(Eigs_plus))
plt.scatter(Eigs_minus,-np.ones_like(Eigs_minus))
