function [eta]=nose_hoover(eta,tau,T,Tset,t_step)

tau2i=1/(tau^2);

eta = eta + t_step*tau2i*((T/Tset) - 1);


