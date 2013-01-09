
%function NanoMD

%Jason Larkin
%September 15, 2009
%Description: MD code to simulate Nanoparticle for MD Class. 
%This code uses units derived from the energy and length scale of the
%Lennard-Jones potential.
%-------------------------------------------------------
%-------------PROGRAM PARAMETERS------------------------
%-------------------------------------------------------
%Declare Global Variables
ncell=10;                            %number of unit cells used in x,y,z
%N=2^(3*ncell-1);                    %number of atoms based on ncells
N=4;                               %number of atoms in simulation
L(1,1:3)=zeros(1,3);                %contains the simulation cell size
pi = atan(1)*4;                            %contains value of pi, in C this would be pi = atan(1) * 4
%a_0 = (3.4E-10)*1.01;
%LJ Potential and Material Parameters
    epsilon_Ar = 1.67E-21;              %aJ (1.67E-21 Joules) aJ=1E-18 J
    sigma_Ar = 3.4E-10;                 %Angstroms 3.4E-10 meters
a_0 = 5.30E-10/sigma_Ar;           %the lattice constant of Ar: http://www.infoplease.com/periodictable.php?id=18
mass_Ar = 6.6326E-26;               %1E-28 kg (6.6326E-26 kg)
mass_Ar = mass_Ar/mass_Ar;
kb = 1.3806E-23;                    %aJ/k (1.3806E-23 J/K)
a = 100;                            %Angstroms cutoff radius in terms of sigma_Ar 
a2 = a*a;
cutoff = 0;                         %true=1, use the cutoff radius. This of course introduces a discontinuity in the force at r=a
%Main Arrays (M X N)
x(1:N,1:3)=zeros(N,3);              %position
p(1:N,1:3)=zeros(N,3);              %momentum
p_half(1:N,1:3)=zeros(N,3);         
x_o(1:N,1:3)=zeros(N,3);            %initial position
m(1,1:N) = zeros(N,1);              %mass of the particle
m(1,1:N) = 1; 
ident_letter(1,1:N)=zeros(N,1);     %letter identifying the particle, e.g. C, N, etc.
rij(1,1:3)=0;                       %pairwise distance between atoms i and j
F(1:N,1:3)=zeros(N,3);           %force on particle (1E-28 kg Anstromg/fs)

%time has units of sigma/sqrt(epsilon/mass)) = secs.  1 time step = 
t_step = 0.002;                    % time step (in units of 2.1427E-12 s)
t_total = 20000;                     % total simulation time (in number of steps at rate t_step)
t = 0;                              % current time (in number of time steps)
t_stats = 10;                       % how often energy and momentum statistics are outputted
props(1,1:7)=0;                     % properties to save
pIndex=1;                           % properties index

%----NOT USED
t_xyz = 1000;                       % how often data is outputed for Chime
t_cfg = t_total+1;                  % how often data is outputed for Atomeye
t_radial = t_total+1;               % how often the radial distribution function is calculated
RDF = zeros(5000,2);
thermostat = 1;	%0=false	
barostat = 0; 			
quench = 0;
eta_t = 0.2;                        % thermostat parameter
tau_t = 0.0;                       % thermostat time constant
%----NOT USED
Tset_K = 20.;                      % temperature set (K) 
Tset = Tset_K*kb/epsilon_Ar;    % reduced temperature
%I don't like this style: * kb * (1./epsilon_Ar); 

%Barostat Params
Pset = 1.;                          % desired simulation pressure here (non-dim)
eps_p = 0.0;                         % barostat parameter
tau_p = 1.0;                        % barostat time constant
V = 0; 
P_viral = 0;
P = 0;
% Quench parameters
eta = 0.;                           % quench parameter, remove kinetic energy from the system
%Calculate Potential Energy cutoff value
pe_cutoff = pe(a2);
%Calculate Force Cutoff value
f_cutoff=force_eval(a2);

%Set initial conditions if the simulation is not a restart
RESTART=0;
%Set to run in NVE to initialize from 0 K state
rescale=0;

%-------------------------------------------------------
%------------INITIALIZATION-----------------------------
%-------------------------------------------------------

%Set initial conditions if the simulation is not a restart
if RESTART == 1    
%[x,v,F,P,P_viral,T,KE,PE]=initialization_fresh(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x,epsilon_Ar,sigma_Ar,m,kb,V,t_step,mass_Ar);
else
    %read_data
    x = dlmread('D:\Classes\CMU\Molecular Simulation\HW2\Q3\10.txt');       %read in initial positions
    x_0=x;                                                  %set initial positions are save
    x_p=x;                                                  %set absolute position if wanting to measure diffusion
    m(1,1:N) = mass_Ar;                                     %set initial masses
    [F,P,P_viral]=force(N,L,a2,f_cutoff,a,x,cutoff,V,0);    %find initial forces
    T=Tset;                                                 %Set temperature as a dummy parameter
    V=1;
    P=1;
    Pset=1;
    L(1,1:3)=1000;                                          %Set simulation domain to basically infinite
end


%-------------------------------------------------------
%------------SIMULATION---------------------------------
%-------------------------------------------------------

for t=1:t_total  
    p_half=half_momentum(N,p,F,t_step,eta_t,eps_p,p_half);
    %[L,V] = baro_size(L,eps_p,t_step,barostat);
    [x,x_p]=position(N,t_step,x,x_p,p_half,p,m,L,tau_p,eps_p,P,Pset);
    %[eps_p,eta_t]=baro_param(barostat,thermostat,eta_t,t_step,tau_t,T,Tset,eps_p,tau_p,P,Pset);
    [F,P,P_viral]=force(N,L,a2,f_cutoff,a,x,cutoff,V,T);
    p=momentum(N,p_half,F,t_step,eta_t,eps_p);
    p_sum=momentum_total(N,p);
    if rem(t,t_stats)==0
        %plot_pos(x,L);
        PE=Pot_Energy(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x);
        KE = ke(N,p,m);
        T=temperature(N,p,m);
        %xrms(pIndex,1)=mean_square(x_p,x_0,N);
        props(pIndex,1) = t;
        props(pIndex,2)=T; props(pIndex,3)=P; props(pIndex,4)=V;
        props(pIndex,5)=KE; props(pIndex,6)=PE; props(pIndex,7)=PE+KE;
        props(pIndex,8)=sum(p(:,1));
        props(pIndex,9)=sum(p(:,2));
        props(pIndex,10)=sum(p(:,3));
        pIndex=pIndex+1;
        %write trajectories to file in CHIME format
        str = 'D:\Classes\CMU\Molecular Simulation\HW2\Q3\props\';
        dlmwrite(strcat(str,'x.xyz'),N-1,'-append');
        dlmwrite(strcat(str,'x.xyz'),x,'-append','delimiter','\t');
    end
end


%-------------------------------------------------------
%------------WRITE PROPERTIES---------------------------
%-------------------------------------------------------

strend='simulation over: ready to write files'
str = 'D:\Classes\CMU\Molecular Simulation\HW2\Q3\props\';
    dlmwrite(strcat(str,'m.dat'),m);
    dlmwrite(strcat(str,'p.dat'),p);
    dlmwrite(strcat(str,'F.dat'),F);
    dlmwrite(strcat(str,'L.dat'),L);
    dlmwrite(strcat(str,'props.dat'),props);


%-------------------------------------------------------
%------------END SIMULATION-----------------------------
%-------------------------------------------------------
    
    


%-------------------------------------------------------
%------------FUNCTIONS----------------------------------
%-------------------------------------------------------

% function phi=pe(r2);
% %FUNCTION: evaluate the potential energy for a given separation r2
%     inv_r6 = 1./(r2*r2*r2);
%     phi = 4.*(inv_r6*inv_r6 - inv_r6);
% end
% 
% 
% function f=force_eval(r2);
% %FUNCTION: evaluate magnitude of the force for given separation r2
%     %using the method derived here:
%     %http://www.pages.drexel.edu/~cfa22/msim/node36.html
%     r4 = r2*r2; r8 = r4*r4; r14i = 1/(r8*r4*r2);
%     f = 48*r14i - 24*r4*r2*r14i;
% end
% 
% function [F,P,P_viral]=force(N,L,a2,f_cutoff,a,x,cutoff,V,T);
% %FUNCTION: find forces acting on all atoms in system.
%     %clear old Forces
%     F(1:N,1:3)=zeros(N,3);
%     P_viral = 0;
%     for i=1:N
%         for j=(i+1):N
%             for k=1:3
%                 rij(1,k) = x(i,k) - x(j,k);
%                 rij(1,k) = periodic(rij(1,k),L(1,k));   %periodic boundaries
%             end
%             r2 = rij(1,1)^2 + rij(1,2)^2 + rij(1,3)^2;
%             if r2<a2
%                 if cutoff == 1
%                     f = force_eval(r2);
%                     %f = force_eval(r2) - f_cutoff*a*(1/sqrt(r2));
%                 else
%                     %f = force_eval(r2) - f_cutoff*a*(1/sqrt(r2));
%                     f = force_eval(r2);
%                 end
%                 for k=1:3    
%                     F(i,k) = F(i,k) + rij(1,k)*f;
%                     F(j,k) = F(j,k) - rij(1,k)*f;
%                 end
%                 P_viral = P_viral +r2*f;
%             end
%         end
%     end
%     %P = (P_viral/(3*V)) + ((N*T)/V);
%     %Set P and Pset to 1 for the purpose of NanoMD
%     P=1;
% end
% 
% function p_half=half_momentum(N,p,F,t_step,eta_t,eps_p,p_half)
% %FUNCTION: find the momentum at the half time step.
%     %p_half(1:N,1:3)=zeros(N,3); 
%     for i=1:N
%         for k=1:3
%             p_half(i,k) = p(i,k) + (F(i,k) - (eta_t+eps_p)*p(i,k))*t_step*0.5;
%         end
%     end
% end
% 
% 
% function [x,x_p]=position(N,t_step,x,x_p,p_half,p,m,L,tau_p,eps_p,P,Pset)
% %FUNCTION: evolve the positions.
%     term1 = eps_p*t_step;
%     term2 = (t_step*t_step/(2*(tau_p*tau_p)))*(P-Pset);
%     term3 = (t_step*eps_p)*(t_step*eps_p)*0.5;
%     for i=1:N
%         for k=1:3
%             x(i,k) = (1 + term1 + term2 + term3)*x(i,k) + p_half(i,k)*t_step/m(1,i) + (t_step*t_step)*eps_p*p(i,k)*0.5/m(1,i);
%             x_p(i,k) = (1 + term1 + term2 + term3)*x_p(i,k) + (t_step*t_step)*eps_p*p(i,k)*0.5/m(1,i);
%             if x(i,k)>L(1,k)   %should it be <=???????
%                 x(i,k) = x(i,k) - L(1,k);
%             end
%             if x(i,k)<0
%                 x(i,k) = x(i,k) + L(1,k);
%             end
%         end
%     end
% end
% 
% function p=momentum(N,p_half,F,t_step,eta_t,eps_p)
% %FUNCTION: calcualte the momentum at the full time step.
%     p(1:N,1:3)=zeros(N,3); 
%     bottom = 1 + (eta_t + eps_p)*t_step*0.5;
%     for i=1:N
%         for k=1:3
%                 p(i,k) = (p_half(i,k) + F(i,k)*0.5*t_step)/bottom;
%         end
%     end
% end
% 
% function p_sum=momentum_total(N,p)
% %FUNCTION: calcualte the total momentum.
%     p_sum(1,1:3)=0;
%     for k=1:3
%         for i=1:N
%             p_sum(1,k) = p_sum(1,k) + p(i,k);
%         end
%     end
% end
% 
% function PE=Pot_Energy(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x)
% %FUNCTION: calcualte the total Potential Energy in the system.
%     PE=0;
%     for i=1:N
%         for j=(i+1):N
%             for k=1:3
%                 rij(1,k) = x(i,k) - x(j,k);
%                 rij(1,k) = periodic(rij(1,k),L(1,k));  %periodic boundaries
%             end
%             r2 = rij(1,1)^2 + rij(1,2)^2 + rij(1,3)^2;
% 
%             if (r2<a2)
%                 if cutoff==1
%                     PE = PE + pe(r2) - pe_cutoff;
%                 else
%                     PE = PE + pe(r2);
%                 end
%             end
%         end
%     end
% end

%end
