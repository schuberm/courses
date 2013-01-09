function [v,T]=berendsen(N,v,m,Tset,kb,tau_t,t_step)

%Using the Berendsen Thermostat which does NOT produce NVT ensemble
%http://www.pages.drexel.edu/~cfa22/msim/node43.html
%I use this for single species.  Multiple species would need v=p./m
mass=m(1,1);
KE=kineticenergy(N,v,m);
T = 2*KE/(3*N*kb);
lambda = sqrt( 1 + (t_step/tau_t)*((T/Tset) - 1));
v = v*lambda;

