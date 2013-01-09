function T=temperature(N,p,m)
KE = ke(N,p,m);
T = (2/3)*((KE)/(N-1));
