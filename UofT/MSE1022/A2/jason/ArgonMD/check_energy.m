function ENERGY=check_energy(t,KE,PE,TE,ENERGY)
%determine periodic boundary conditions

ENERGY(t,1)=KE;
ENERGy(t,2)=PE;
ENERGY(t,3)=TE;

