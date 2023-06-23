N=150; L=40; Nx=4*N; dx=L/Nx; x=-L/2:dx:L/2-dx; k0=2*pi/L;
mu=1; B=-1/sqrt(1+16*mu/3);
u=sqrt(-4*B*mu./(B+cosh(2*sqrt(mu)*x)));

G1=-mu-u.^2+u.^4; 
G2=-mu-3*u.^2+5*u.^4;

for n=-N:1:N
c1(n+N+1)=dx*sum(G1.*exp(-i*k0*n*x))/L;
c2(n+N+1)=dx*sum(G2.*exp(-i*k0*n*x))/L;
end

D=(i*k0*diag([-N:N])).^2;
C0=zeros(2*N+1);
C1=toeplitz([c1(N+1:2*N+1) zeros(1,N)],[c1(N+1:-1:1) zeros(1,N)]);
C2=toeplitz([c2(N+1:2*N+1) zeros(1,N)],[c2(N+1:-1:1) zeros(1,N)]);
M=[ C0
D+C1
D+C2
-C0 ];
eigvalues=eig(i*M);
plot(eigvalues, '.'); axis([-5 5 -5 5]);
