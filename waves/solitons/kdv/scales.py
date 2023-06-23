


def rescale(x,u,w,r_out):

    r = 7

    phi3 = -1
    phi5 = -2.1

    wt = w/2

    chi = r_out/r

    bond = 0.85**2/wt**2/chi**4

    delta = abs(-phi3/6-bond)

    gamma = (-phi5/120-bond*phi3/6+delta**2/4)

    epsilon2 = 2*chi**2*gamma/delta**2

    x_scale = wt*chi/r_out*np.sqrt(delta/2)
    a_scale = 4*wt/5


    return u*a_scale,x*x_scale


