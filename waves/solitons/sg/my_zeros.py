
def M12_finder(E):
    
    Mt = B(u[n],ux,ut,E+0j)
    
    return Mt[0,1].real

def trace_finder_plus(E):
    Mt = B(u[n],ux,ut,E+0j)

    tr = np.trace(Mt.real)/2

    if np.abs(tr)<=1:
        return tr-1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))-1
    
def trace_finder_minus(E):
    Mt = B(u[n],ux,ut,E+0j)

    tr = np.trace(Mt.real)/2
    
    if np.abs(tr)<=1:
        return tr+1

    else:
        return np.sign(tr)*(1+np.log(np.abs(tr)))+1


Eigs_plus = []
Eigs_minus = []

mins = argrelmin(trace)[0]
maxs = np.concatenate(([int(mins[0]/2)],argrelmax(trace)[0]))


Emaxs = E[maxs]
Emins = E[mins]



Ebr = np.zeros(int(len(Emaxs)+len(Emins)))

for i in range(len(Emaxs)):

    Ebr[2*i] = Emaxs[i].real

for i in range(len(Emins)):
    
    Ebr[2*i+1] = Emins[i].real



for i in range(len(Ebr)-1):

    if (np.sign(trace_finder_plus(Ebr[i]))!=np.sign(trace_finder_plus(Ebr[i+1]))):
    
        Eigs_plus.append(bisect(trace_finder_plus,Ebr[i],Ebr[i+1],xtol=1e-16))

    if (np.sign(trace_finder_minus(Ebr[i]))!=np.sign(trace_finder_minus(Ebr[i+1]))):
        
        Eigs_minus.append(bisect(trace_finder_minus,Ebr[i],Ebr[i+1],xtol=1e-16))


#Eref =  Eigs_minus[2]

#Eigs_plus = np.array(Eigs_plus)-Eref
#Eigs_minus = np.array(Eigs_minus)-Eref
#Er = E-Eref

#Ems = E[mins]
#Ems = Ems[Ems>0]

#Eigs_plus = np.concatenate((Eigs_plus,E[maxs]))
#Eigs_minus = np.concatenate((Eigs_minus,Ems))


plt.figure()
plt.title("Trace")
plt.plot(E,trace)
#plt.plot(E,S)
plt.axhline(1,0,1,ls="--",color="tab:red")
plt.axhline(-1,0,1,ls="--",color="tab:red")
plt.axvline(0,0,1,ls="--",color="tab:red")
plt.axhline(0,0,1,ls="--",color="tab:red")

plt.xlabel("E")
plt.ylabel("Trace")

plt.scatter(Eigs_plus,np.ones_like(Eigs_plus))
plt.scatter(Eigs_minus,-np.ones_like(Eigs_minus))


Eigs = np.concatenate((Eigs_plus,Eigs_minus))
Eigs = np.sort(Eigs)

sindex = []

for i in range(1,int(len(Eigs)/2)):
    sindex.append((Eigs[2*i]-Eigs[2*i-1])/(Eigs[2*i]-Eigs[2*i-2]))


j = np.arange(len(sindex))
plt.figure()
plt.scatter(j,sindex)


def amps(jr):
    
    Eref = Eigs[2*jr]
    amps = 2*(Eref-Eigs_minus)/lam

    return amps
