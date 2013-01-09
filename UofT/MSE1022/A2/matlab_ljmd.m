function sys = matlab_ljmd()

% 
% simple lennard-jones potential MD code with velocity verlet.
% units: Length=Angstrom, Mass=amu; Energy=kcal
% 
% Code ported from Axel Kohlmeyer's Parallel LJMD code availiable 
% at http://sites.google.com/site/akohlmey/software/ljmd
% 
% I've retained as much of the variable naming conventions and
% format as possible so as to be comparable to the original code
% and easily updated / changed. However, in changing to matlab
% there are a number of optimized vector based functions that could likely
% improve the performace. I leave that to the user.
% -Chris MacDermaid

% To run, just make sure you are in the directory containg this
% file, and the restart file. Make sure you have the correct box
% dimensions and number of atoms and type: matlab_ljmd


 %% Change format to long double precision
 format long;

 %% Declare some necessary constants and globals
 global kboltz mvsq2e sys
 
 kboltz = 0.0019872067; % Boltzman constant in kcal/mol/K
 mvsq2e = 2390.05736153349; % m*v^2 in kcal/mol  
 
 %% Declare some variables for tracking energetics
 sys.ekin = 0.0;                    % Kinetic Energy
 sys.epot = 0.0;                    % Potential Energy
 sys.temp = 0.0;                    % System Temperature (microcanonical ensemble)
 sys.nfi = 1;                       %Step Counter
 
 %%%%%%%%%%%%%%%%%%%%%
 %% Begin User Input%%
 %%%%%%%%%%%%%%%%%%%%%
 
 %% System Dependent
 sys.natoms = 108;                  % Number of atoms
 sys.mass = 39.948;                 % Argon Mass in AMU
 sys.epsilon = 0.2379;              % epsilon in kcal/mol
 sys.sigma = 3.405;                 % sigma in angstrom
 sys.rcut = 12.0;                   % rcut in angstrom
 sys.box = 17.1580;                 % Box Length in angstrom
 sys.nsteps = 10000;                % number of MD timesteps
 sys.dt = 5.0;                      % Timesteps in fs
 nprint = 100;                      % Output Frequency (number of MD steps)
 
 %% Output Files
 sys.restfile = 'argon_108.rest'; %Restart File
 sys.trajfile = 'argon_108.xyz';  %Trajectory
 sys.ergfile = 'argon_108.dat';   %energies

 %%%%%%%%%%%%%%%%%%%
 %% End User Input%%
 %%%%%%%%%%%%%%%%%%%
 
 %% Initilize arrays for storage.
 
 % Position
 sys.rx = zeros(1,sys.natoms);
 sys.ry = zeros(1,sys.natoms);
 sys.rz = zeros(1,sys.natoms);
 
 % Velocity
 sys.vx = zeros(1,sys.natoms);
 sys.vy = zeros(1,sys.natoms);
 sys.vz = zeros(1,sys.natoms);
 
 % Force
 sys.fx = zeros(1,sys.natoms);
 sys.fy = zeros(1,sys.natoms);
 sys.fz = zeros(1,sys.natoms);
 
 %% Read the restart file:
 fid=fopen(sys.restfile, 'r');
 [A,count] = fscanf(fid,'%f %f %f', [3 inf]);
 fclose(fid);
 
%% Fill the r and v arrays with the restart info
 sys.rx = A(1,[1:sys.natoms]);
 sys.ry = A(2,[1:sys.natoms]);
 sys.rz = A(3,[1:sys.natoms]);
 sys.vx = A(1,[sys.natoms + 1:sys.natoms * 2]);
 sys.vy = A(2,[sys.natoms + 1:sys.natoms * 2]);
 sys.vz = A(3,[sys.natoms + 1:sys.natoms * 2]);
 

 %% Open the log and trajectory files for writing
 erg = fopen(sys.ergfile, 'w');
 traj = fopen(sys.trajfile, 'w');

 %% Initilize our values of force / kinetic energy and output
 %initial values
 force();
 ekin();

 fprintf('Starting simulation with %d atoms for %d steps.\n',sys.natoms, sys.nsteps);
 fprintf('     NFI            TEMP            EKIN                 EPOT              ETOT\n');
 output(erg,traj);
 
 %% The Main MD Loop
 for i = 1 : 1 : sys.nsteps
   sys.nfi = i;
   
   %% Output? 
   if mod(i,nprint) == 0; output(erg,traj); end; 
   
   %Propogate and recompute Energies
   velverlet();
   ekin();
      
 end  
 
 %% Cleanup
 fclose(erg);
 fclose(traj);
 
 % Minimum image convention helper
 function x = pbc(x, boxby2, box)
   while x > boxby2; x = x - box ;end 
   while x < -boxby2; x = x + box;end
 
 % Compute the kinetic energy of our system  
 function ekin()
  global kboltz mvsq2e sys

  %% Clear kinetic energy
  sys.ekin = 0;
  
  for i = 1 : 1 : sys.natoms
    sys.ekin = sys.ekin + sys.vx(i)*sys.vx(i) ...
        + sys.vy(i)*sys.vy(i) ...
        + sys.vz(i)*sys.vz(i);
  end
  
  sys.ekin = sys.ekin * 0.5 * mvsq2e * sys.mass;
  sys.temp = 2.0 * sys.ekin / (3.0 * sys.natoms - 3.0) / kboltz;
  

% Compute Forces on Atoms  
function force()

 global sys 
 
 %% Define some variables for the force calculation
 %% ahead of time to save on computations
 ffac = 0.0;
 c12 = 4.0 * sys.epsilon * sys.sigma ^ 12;
 c6 = 4.0 * sys.epsilon * sys.sigma ^ 6;
 boxby2 = 0.5 * sys.box;
 rcsq = sys.rcut * sys.rcut;
 
 %% zero out our force arrays
 sys.fx = zeros(1,sys.natoms);
 sys.fy = zeros(1,sys.natoms);
 sys.fz = zeros(1,sys.natoms);
 
 %% Clear the potential energy
 sys.epot = 0.0;
 
 for i = 1 : 1 : sys.natoms - 1
   rx1 = sys.rx(i);   
   ry1 = sys.ry(i);
   rz1 = sys.rz(i);
   
   for j = i + 1 : 1 : sys.natoms

     rx = rx1 - sys.rx(j);
     ry = ry1 - sys.ry(j);
     rz = rz1 - sys.rz(j);

     % apply PBC
     while rx >  boxby2; rx = rx - sys.box; end 
     while rx < -boxby2; rx = rx + sys.box; end
     while ry >  boxby2; ry = ry - sys.box; end 
     while ry < -boxby2; ry = ry + sys.box; end
     while rz >  boxby2; rz = rz - sys.box; end 
     while rz < -boxby2; rz = rz + sys.box; end
     
     rsq = rx * rx + ry * ry + rz * rz;
 
     if rsq < rcsq
       rinv = 1.0 / rsq;
       r6 = rinv * rinv * rinv;
       
       ffac = (12.0 * c12 * r6 - 6.0 * c6) * r6 * rinv;
       sys.epot = sys.epot + r6 * (c12 * r6 - c6);
       
       %% Apply Newton's third law
       sys.fx(i) = sys.fx(i) + rx*ffac;
       sys.fy(i) = sys.fy(i) + ry*ffac; 
       sys.fz(i) = sys.fz(i) + rz*ffac; 
       sys.fx(j) = sys.fx(j) - rx*ffac; 
       sys.fy(j) = sys.fy(j) - ry*ffac; 
       sys.fz(j) = sys.fz(j) - rz*ffac; 
     end
   
   end
 
 end  
 
%% Velocity Verlet
function velverlet()
   
   global mvsq2e sys

   dtmf = 0.5 * sys.dt / mvsq2e / sys.mass;

   %% Propagate velocities by half step and positions by full step
   for i = 1 : 1 : sys.natoms
     % Vels
     sys.vx(i) = sys.vx(i) + dtmf * sys.fx(i);
     sys.vy(i) = sys.vy(i) + dtmf * sys.fy(i);
     sys.vz(i) = sys.vz(i) + dtmf * sys.fz(i);
     
     %Pos
     sys.rx(i) = sys.rx(i) + sys.dt * sys.vx(i);
     sys.ry(i) = sys.ry(i) + sys.dt * sys.vy(i);
     sys.rz(i) = sys.rz(i) + sys.dt * sys.vz(i);
   end
   
   %% Compute our forces and PE
   force();
      
   %% Propagate velocities for the remainining half step
   for i = 1 : 1 : sys.natoms
     sys.vx(i) = sys.vx(i) + dtmf * sys.fx(i);
     sys.vy(i) = sys.vy(i) + dtmf * sys.fy(i);
     sys.vz(i) = sys.vz(i) + dtmf * sys.fz(i);
    end  


%% append data to output
function output(erg, traj)
      global sys  

      fprintf(1,'% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n', ...
                  sys.nfi, sys.temp, sys.ekin, sys.epot, sys.ekin + ...
                  sys.epot)
      
      fprintf(erg, '% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n', sys.nfi, ...
                    sys.temp, sys.ekin, sys.epot, sys.ekin + sys.epot);
      
      fprintf(traj, '%d\n nfi=%d etot=%20.8f\n', ...
              sys.natoms, sys.nfi, sys.ekin + sys.epot);
      
      for i = 1 : 1 : sys.natoms 
        fprintf(traj, 'Ar  %20.8f %20.8f %20.8f\n', sys.rx(i), ...
                sys.ry(i), sys.rz(i));
      end  
     
 % Plots
 figure(1);
 data=load('argon_108.dat');
 subplot(2,2,1); plot(data(:,1),data(:,2)); grid on; xlabel('Timestep');ylabel('Temperature');
 subplot(2,2,2); plot(data(:,1),data(:,3)); grid on; xlabel('Timestep');ylabel('Kinetic energy');
 subplot(2,2,3); plot(data(:,1),data(:,4)); grid on; xlabel('Timestep');ylabel('Potential energy');
 subplot(2,2,4); plot(data(:,1),data(:,5)); grid on; xlabel('Timestep');ylabel('Total energy');
