function [x,v,F,P,P_viral,T,KE,PE,V]=initialization_fresh(thermostat,barostat,N,L,pe_cutoff,a,a2,f_cutoff,cutoff,m,t_step,mass_Ar,a_0,ncell,x,F,v,eta_t,eps_p,Tset,Pset)

%-------------------------------------------------------
%------------SIMULATION INITIALIZATION------------------
%-------------------------------------------------------

%Produce Intial lattice
out = fcclattice(ncell,ncell,ncell); x = a_0*(out-1);
%Set initial species mass
m(1,1:N) = mass_Ar; 
%Randomize Initial Coordiantes away from Equillirbium by 5% a_0
% latper = 0.05;
% x = x +a_0*latper*(rand(N,3)-0.5);
% %Adjust the positions so that all random coords are inside box
% x = x + abs(min(min(x)))+a_0*latper;
% x_0=x;
% simulation cell size
	for k=1:3 
        L(1,k) = max(max(x))+a_0/2;	
        % desired system size. May be specified here, determined from input coords.
        %So, I add an extra a_0 so that using the periodic boundary
        %conditions adds an extra lattice constant in between mirrored
        %particles.
    end
	V = L(1,1)*L(1,2)*L(1,3);    %system volume
    plot_pos(x,L);

            PE=Pot_Energy(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x);
            KE = ke(N,v,m);
            T=temperature(N,v,m);
            % find initial forces/pressure
            [F,P,P_viral]=initial_force(N,F,L,a2,f_cutoff,a,x,cutoff,V,T);
            % find initial barostat and thermostat parameters
            %baro_initial();
            %Set Initial velocities Randomly
            v = initial_vel(N,v);
            v=vel_rescale(N,v,m,Tset);
            
            %Thermostat Parameters
            if thermostat == 1 
                eta_t = eta_t + t_step*(1./(tau_t*tau_t))*((T/Tset)-1.);
            else
                eta_t = 0.;
            end
            
            %Barostat Parameters
            if barostat == 1 
                eps_p = eps_p + t_step*(1./(tau_p*tau_p))*(P-Pset);
            else
                eps_p = 0.;
            end
            
            
