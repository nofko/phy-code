


ro = np.linspace(8,8.36,500)
ro2 = np.linspace(8.52,8.9,500)
r = 7

phi3 = -1
phi5 = -2.1

def epsilon2(r_out,w):

    #w = 2*(r_out-r)
    wt = w/2

    chi = r_out/r
    
    bond = 0.85**2/wt**2/chi**4

    delta = -phi3/6-bond

    gamma = (-phi5/120-bond*phi3/6+delta**2/4)
    
    result = 2*chi**2*gamma/delta**2

    print(delta)
    
    return result



plt.figure(figsize=[12,9])
#plt.plot(w,k3,label="$k^4$")
#plt.plot(w,k5,label="$k^6$")

plt.plot(ro,epsilon2(ro))
plt.plot(ro2,epsilon2(ro2))

plt.axhline(0,linestyle="--")

plt.xlabel("$R_o$ [cm]",fontsize=35)
plt.ylabel("$\epsilon^2$",fontsize=35)
#plt.legend()
