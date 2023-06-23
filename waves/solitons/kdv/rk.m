function [tdata udata] = KdVsolu(uex, N, tmax, dt, nplot)
x = (2*pi/N)*(-N/2:N/2-1)’;
u = feval(uex, x, 0); v = fft(u);
k = [0:N/2-1 0 -N/2+1:-1]’; ik3 = 1i*k.^3;
nplt = floor((tmax/nplot)/dt); nmax = round(tmax/dt);
udata = u; tdata = 0;
for n = 1:nmax
t = n*dt; g = -.5i*dt*k;
E = exp(dt*ik3/2); E2 = E.^2;
a = g.*fft(real( ifft( v ) ).^2);
b = g.*fft(real( ifft(E.*(v+a/2)) ).^2); % 4th-order
c = g.*fft(real( ifft(E.*v + b/2) ).^2); % Runge-Kutta
d = g.*fft(real( ifft(E2.*v+E.*c) ).^2);
v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;
if mod(n,nplt) == 0 u = real(ifft(v));
udata = [udata u]; tdata = [tdata t];
end
end
end