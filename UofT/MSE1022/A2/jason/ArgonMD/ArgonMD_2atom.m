
%Jason Larkin
%September 15, 2009
%Description: MD code to simulate Argon. 
%Code stolen from Eric Landry.
%This code uses all normal metric units, e.g. meters, seconds, kg Joules,
%etc.
%Future versions will used some sort of modified units.
%-------------------------------------------------------
%-------------PROGRAM PARAMETERS------------------------
%-------------------------------------------------------
%Declare Global Variables
ncell=2;                            %number of unit cells used in x,y,z
%N=2^(3*ncell-1);                    %number of atoms based on ncells
N=2;
L(1,1:3)=zeros(1,3);                %contains the simulation cell size
pi = pi;                            %contains value of pi, in C this would be pi = atan(1) * 4
%a_0 = 5.260E-10;                    %the lattice constant of Ar: http://www.infoplease.com/periodictable.php?id=18
%a_0 = (3.4E-10)*1.01;
a_0 = 10E-11;
%Specify LJ values for Argon
epsilon_Ar = 1.67E-21;              %aJ (1.67E-21 Joules) aJ=1E-18 J
sigma_Ar = 3.4E-10;                 %Angstroms 3.4E-10 meters
a_0=(2^(1/6))*sigma_Ar;
mass_Ar = 6.6326E-26;               %1E-28 kg (6.6326E-26 kg)
kb = 1.3806E-23;                    %aJ/k (1.3806E-23 J/K)
a = 2.5*sigma_Ar;                   %Angstroms cutoff radius in terms of sigma_Ar 
a2 = a*a;
cutoff = 1;                         %true=1, use the cutoff radius. This of course introduces a discontinuity in the force at r=a
%Main Loop Control Variable
x(1:N,1:3)=zeros(N,3);              %position
v(1:N,1:3)=zeros(N,3);              %velocity
x_o(1:N,1:3)=zeros(N,3);            %initial position
m(1,1:N) = zeros(N,1);              %mass of the particle
ident_letter(1,1:N)=zeros(N,1);     %letter identifying the particle, e.g. C, N, etc.
rij(1,1:3)=0;                       %pairwise distance between atoms i and j
F(1:N,1:3)=zeros(1:N,1:3);           %force on particle (1E-28 kg Anstromg/fs)

t_step = 0.1E-15;                     % time step (in s)
t_total = 5000;                     % total simulation time (in number of steps at rate t_step)
t = 0;                              % current time (in number of time steps)
t_stats = 10;                       % how often energy and momentum statistics are outputted

%----NOT USED
t_xyz = 1000;                       % how often data is outputed for Chime
t_cfg = t_total+1;                  % how often data is outputed for Atomeye
t_radial = t_total+1;               % how often the radial distribution function is calculated
thermostat = 0;	%0=false	
barostat = 0; 			
quench = 0;
eta_t = 0.;                         % thermostat parameter
tau_t = 0.05;                       % thermostat time constant
%----NOT USED
Tset_K = 40.;                       % temperature set (K) 
Tset = Tset_K;                      % temperature to set system to (K) 
%I don't like this style: * kb * (1./epsilon_Ar); 
Pset = 0.;                          % desired simulation pressure here (non-dim)
eps_p = 0.;                         % barostat parameter
tau_p = 1.0;                        % barostat time constant
% Quench parameters
eta = 0.;                           % quench parameter, remove kinetic energy from the system

%-------------------------------------------------------
%------------SIMULATION INITIALIZATION------------------
%-------------------------------------------------------
%Produce 2 Argon atoms and check force without cutoff
% x(1,1:3)=0;
% a1 = 3E-10; af = 5E-10; cnt=1;
% for aplot=a1:(af-a1)/1000:af
%     delta = a/sqrt(3);
%     x(2,1)=aplot;
%     x(2,2)=0; x(2,3)=0;
%     for k=1:3
%             rij(1,k) = x(1,k) - x(2,k);
%             %rij(1,k) = periodic(rij(1,k),L(1,k));   %periodic boundaries
%         end
%     r2 = rij(1,1)^2 + rij(1,2)^2 + rij(1,3)^2;
%     PE(cnt,1)=pe(r2,epsilon_Ar,sigma_Ar);
%     f = force_eval(r2,epsilon_Ar,sigma_Ar);
%     F(cnt,1) = rij(1,1)*f;
%     F(cnt,2) = rij(1,2)*f;
%     F(cnt,3) = rij(1,3)*f;
%     cnt=cnt+1;
% end
% aplot=a1:(af-a1)/1000:af;
% plot(aplot,PE)
% plot(aplot,F(:,1))
% pause

%Check Force with cutoff on 2 Argon Atoms
f_cutoff=force_eval(a2,epsilon_Ar,sigma_Ar);
V=1; %set V to dummy value of 1
x(1,1:3)=0;
a1 = 3E-10; af = 10E-10; cnt=1;
for aplot=a1:(af-a1)/1000:af
    delta = a/sqrt(3);
    x(2,1)=aplot;
    x(2,2)=0; x(2,3)=0;
    [F,P,P_viral,v]=force(N,F,L,a2,f_cutoff,a,x,cutoff,epsilon_Ar,sigma_Ar,v,t_step,V,Tset);
    Fcalc(cnt,1) = F(1,1);
    Fcalc(cnt,2) = F(1,2);
    Fcalc(cnt,3) = F(1,2);
    cnt=cnt+1;
end
aplot=a1:(af-a1)/1000:af;
plot(aplot,Fcalc(:,1))
pause











%Produce Intial lattice
out = fcclattice(ncell,ncell,ncell); x = 2*a_0*(out-min(min(out)));
%Set initial species mass
m(1,1:N) = mass_Ar; 
%Randomize Initial Coordiantes away from Equillirbium by 5% a_0
latper = 0.05;
x = x +a_0*latper*(rand(N,3)-0.5);
%Adjust the positions so that all random coords are inside box
x = x + abs(min(min(x)))+a_0*latper;
%Set Initial velocities Randomly
v = initial_vel(N,v,Tset,kb,mass_Ar);
% simulation cell size
	for k=1:3 
        L(1,k) = max(max(x))+a_0;	
        % desired system size. May be specified here, determined from input coords.
        %So, I add an extra a_0 so that using the periodic boundary
        %conditions adds an extra lattice constant in between mirrored
        %particles.
    end
	V = L(1,1)*L(1,2)*L(1,3);    %system volume
    %Calculate Potential Energy cutoff value
    pe_cutoff = pe(a2,epsilon_Ar,sigma_Ar);
    %Calculate Force Cutoff value
    f_cutoff=force_eval(a2,epsilon_Ar,sigma_Ar);
pause

RESTART='false';
if RESTART == 'false'    
        PE=Pot_Energy(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x,epsilon_Ar,sigma_Ar);
        KE = kineticenergy(N,v,m);
        T=temperature(KE,N,kb);
		% find initial forces/pressure
        [F,P,P_viral,v]=force(N,F,L,a2,f_cutoff,a,x,cutoff,epsilon_Ar,sigma_Ar,v,t_step,V,T);
		% find initial barostat and thermostat parameters
        %baro_initial();
end
   
%-------------------------------------------------------
%------------EQUILIBRIATION-----------------------------
%-------------------------------------------------------
		for t=1:t_total
            [x]=position(N,barostat,t_step,x,F,v,m,L);
            % find forces/pressure/velocity at t+delt
            [F,P,P_viral,v]=force(N,F,L,a2,f_cutoff,a,x,cutoff,epsilon_Ar,sigma_Ar,v,t_step,V,T);
            %Use simple velocity rescaling
            v=vel_rescale(N,v,m,Tset,kb);
            if rem(t,t_stats)==0 
                plot_pos(x,L);
            pause
            PE=Pot_Energy(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x,epsilon_Ar,sigma_Ar)
            KE = kineticenergy(N,v,m)
            T=temperature(KE,N,kb)
            pause
            end
        end

%-------------------------------------------------------
%------------SIMULATION---------------------------------
%-------------------------------------------------------
eIndex=1;
		for t=1:t_total
			%One should have F(t) from the
            %initalization above.  Remember, if one wants to perform
            %restarts then the forces F(t) must be saved to be ablve to
            %immediately jump back into this loop.
            %output to screen every t_stats
            if rem(t,t_stats)==0 
            plot3(x(:,1),x(:,2),x(:,3),'.')
            %pause
            end
            KE = kineticenergy(N,v,m);
            %PROGRESS: what about factor of 2 T = (2/3)*((2*KE)/(N*kb)
            % scale momentum to get close to desired temperature
            [v,T]=berendsen(N,v,m,Tset,kb,tau_t,t_step);
            %find positions at t+t_step. 
            [x]=position(N,barostat,t_step,x,F,v,m,L);
            % find forces/pressure/velocity at t+delt
            [F,P,P_viral,v]=force(N,F,L,a2,f_cutoff,a,x,cutoff,epsilon_Ar,sigma_Ar,v,t_step,V,T);
            %find velocities 
            %[v]=velocity(F,F_half,x,t_step,m);
            %Apply barostat to find new cell size, if active and calculate
            %new barostat parameters
            %[L,V,eps_p,eta_t]=baro_size(N,barostat,L,eps_p,t_step,V,eta_t,tau_t,T,Tset,P,Pset,thermostat);	  
			% find new temperature
            KE=kineticenergy(N,v,m);
			% find new pressure
			P = (P_viral/(3*V)) + ((N*T)/V);
% Calculate PE, KE, TE and display				
			if rem(t,t_stats)==0
            [KE,PE,TE]=total_energy(N,x,L,cutoff,f_cutoff,a,pe_cutoff,v,m,a2,epsilon_Ar,sigma_Ar);
                ENERGY(eIndex,1)=KE;
                ENERGY(eIndex,2)=PE;
                ENERGY(eIndex,3)=TE;
                eIndex=eIndex+1;
            end
									
        end






