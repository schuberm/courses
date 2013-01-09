
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
ncell=4;                            %number of unit cells used in x,y,z
%N=2^(3*ncell-1);                    %number of atoms based on ncells
N=256;
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
%Calculate Potential Energy cutoff value
pe_cutoff = pe(a2,epsilon_Ar,sigma_Ar);
%Calculate Force Cutoff value
f_cutoff=force_eval(a2,epsilon_Ar,sigma_Ar);
a1 = 2E-10; a2 = 10E-10; cnt=1; NUM=100;
for a=a1:(a2-a1)/NUM:a2
    %Produce Intial lattice
    out = fcclattice(ncell,ncell,ncell); x = a*(out-min(min(out)));
    %Set initial species mass
    m(1,1:N) = mass_Ar; 
    % simulation cell size
        for k=1:3 
            L(1,k) = max(max(x))+a_0;	
            % desired system size. May be specified here, determined from input coords.
            %So, I add an extra a_0 so that using the periodic boundary
            %conditions adds an extra lattice constant in between mirrored
            %particles.
        end
        V = L(1,1)*L(1,2)*L(1,3);    %system volume
        PE(cnt,1)=Pot_Energy(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x,epsilon_Ar,sigma_Ar)
        cnt=cnt+1;
        
end
a=a1:(a2-a1)/NUM:a2;
plot(a,PE)
pause




