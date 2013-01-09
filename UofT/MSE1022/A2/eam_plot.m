% Matlab/octave script to generate potential energy vs. atomic separation plot for EAM potential
clear; clc;
% Data from EAM potential file
Nr=500; Nrho=500; dr=0.0112; rcut=5.5500000000000114e+00;
eam=load('pot.data');

F=eam(1:100,:); Frho=reshape(F',500,1);
rp=eam(101:200,:); rphi=reshape(rp',500,1);
rho1=eam(201:300,:); rho=reshape(rho1',500,1);

Nr1=450;rho1=rho(1:Nr1);Frho1=Frho(1:Nr1);
for i=1:Nr1
    r(i)=dr*i;
    phi(i)=rphi(i)/r(i);
    Fr(i)=interp1(rho1,Frho1,rho(i));
end
pot=Fr+0.5*phi;

ep=0.458;
s=2.569;

x=0:0.01:10;
v=4.*ep.*((s./x).^12-(s./x).^6);
f=-4*ep.*s.^6*6./x.^7+4*ep.*s.^12*6.*12./x.^13;
% Plots

figure(2); % zoomed in figures
plot(r,pot,x,v,'--','LineWidth',2); xlim([0.1 5]); ylim([-20 50]); grid on; xlabel('atomic separation (Angstrom)'); ylabel('EAM potential, phi (eV)');legend('EAM','LJ');
print('-dpsc','fig2.eps');

% Force plot (Very approximate)
force=-diff(pot)./diff(r);
figure(3);
plot(r,pot,x,f,'--','LineWidth',2); xlim([0.1 5]);ylim([-20 50]); grid on; xlabel('atomic separation (Angstrom)'); ylabel('Interatomic force (eV/A)');legend('EAM','LJ');
print('-dpsc','fig3.eps');

% % Fit a polynomial curve to the potential energy curve
% Note: poly fit doesnt work well
% byr=1./r;
% [p,S,mu] = polyfit(byr,pot,5); %p=p1*r5+p2*r4+p3*r3+p4*r2+p5*r+p6;
% p_val=polyval(p,byr);
% figure(3); plot(byr,pot,byr,p_val,'--');
% figure(4); plot(r,pot,r,p_val,'--'); xlim([0.1 1]);


