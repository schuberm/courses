function [x,v,F,P,P_viral,T,KE,PE]=initialization(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x,m,V,t_step,mass_Ar,v)

%-------------------------------------------------------
%------------SIMULATION INITIALIZATION------------------
%-------------------------------------------------------

%Set initial species mass
m(1,1:N) = mass_Ar; 
%Find the system length parameters and volume
	for k=1:3 
        L(1,k) = max(max(x))+a_0/2;	
        % desired system size. May be specified here, determined from input coords.
        %So, I add an extra a_0 so that using the periodic boundary
        %conditions adds an extra lattice constant in between mirrored
        %particles.
    end
	V = L(1,1)*L(1,2)*L(1,3);    %system volume
    %calculate initial energies and termperature
    PE=Pot_Energy(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x);
    KE = ke(N,v,m);
    T=temperature(N,v,m);
    [F,P,P_viral,v]=force(N,L,a2,f_cutoff,a,x,cutoff,v,V,T,mass);
    
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
