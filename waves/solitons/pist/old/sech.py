
Eigs_plus = []
Eigs_minus = []


for i in range(len(mus)-1):

    if (np.sign(trace_finder_plus(mus[i]))!=np.sign(trace_finder_plus(mus[i+1]))):

        Eigs_plus.append(bisect(trace_finder_plus,mus[i],mus[i+1],xtol=1e-12))

    if (np.sign(trace_finder_minus(mus[i]))!=np.sign(trace_finder_minus(mus[i+1]))):
        
        Eigs_minus.append(bisect(trace_finder_minus,mus[i],mus[i+1],xtol=1e-12))



maxs = argrelmax(trace)[0]
mins = argrelmin(trace)[0]

Eigs_plus = np.array(Eigs_plus)
Eigs_minus = np.array(Eigs_minus)

Ems = E[mins]
Ems = Ems[Ems>0]

Eigs_plus = np.concatenate((Eigs_plus,E[maxs]))
Eigs_minus = np.concatenate((Eigs_minus,Ems))


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


Eigs = np.concatenate((Eigs_plus,Eigs_minus))
Eigs = np.sort(Eigs)

sindex = []

for i in range(1,int(len(Eigs)/2)):
    sindex.append((Eigs[2*i]-Eigs[2*i-1])/(Eigs[2*i]-Eigs[2*i-2]))


j = np.arange(len(sindex))
    
plt.figure()
plt.scatter(2*np.pi*j/L,sindex)
