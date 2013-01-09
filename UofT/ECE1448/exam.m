
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
	c=299792458;
	zo=(4*pi*10^-7./8.854187817620e-12)^0.5
	eV=6.24150974E18
for i = 1:l
    for j = 1:l
	for k = 1:l
		W(i,j,k)= i^2*pi^2*hbar^2/(2*me*b^2)+j^2*pi^2*hbar^2/(2*me*b^2)+k^2*pi^2*hbar^2/(2*me*b^2);
		E=E+exp(-W(i,j,k)./(kb*T));	
	end
    end
end
    E
    W(1,1,1)*eV
    W(1,1,2)*eV
    (W(1,1,2)-W(1,1,1))/hbar
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
    n(1,1,2)
%Q1
   matele=(8./(5^3)*-200.*5./2*5./2./(9*pi^2))./2^0.5
   matele1=(8./(b^3)*1/400000000*1/400000000*-1/(45000000000000000*pi^2))
   matelel2=-80/(9*pi^2)./2^0.5
   omega=5000E6;
   % Pabs from notes:	
   Ratio=0.25*zo/hbar^2*q^2*matele1^2.*1./(omega)*(n(1,1,1)-3*n(1,1,2))*hbar*(W(1,1,2)-W(1,1,1))./hbar/100
   1./omega
   hbar*(W(1,1,2)-W(1,1,1))./hbar
   
%*(b)*(2)^0.5
	
%   I=c./(8*pi)
   alpha=10*log10(e)*Ratio
%.*(b*100)*(2)^0.5)

   	
%Q3
   mu=9.27400915E-24;
   bohrmag=q*hbar/(2*me)
   omega2=2*0.1*mu/hbar/(2*pi)
   q*0.1/me/(2*pi)
