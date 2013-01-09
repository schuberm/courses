
	clear all
	close all
 	clc
	kb=1.3806488E-23;
	T=300;
	l=5;
	hbar=1.054571726E-34;
	b=5E-9;
	me=9.10938188E-31;
	W=zeros(l,l,l);
	E=0;
	q=1.602E-19;
	c=29979245800;
	zo=(4*pi*10^-7./8.854187817620e-12)^0.5
	eV=6.24150974E18
for i = 1:2
    for j = 1:2
	for k = 1:2
		W(i,j,k)= i^2*pi^2*hbar^2/(2*me*b^2)+j^2*pi^2*hbar^2/(2*me*b^2)+k^2*pi^2*hbar^2/(2*me*b^2);
		E=E+exp(-W(i,j,k)./(kb*T));	
	end
    end
end
    N=1.E16;
    n=zeros(l,l,l);
for i = 1:l
    for j = 1:l
	for k = 1:l
		n(i,j,k)= N*exp(-W(i,j,k)/(kb*T))./E;	
	end
    end
end
    n(1,1,1)
    3*n(1,1,2)

