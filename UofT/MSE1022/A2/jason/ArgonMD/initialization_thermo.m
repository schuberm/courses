function [eta_t,eps_p]=initialization_thermo(thermostat,barostat,t_step,eta_t,eps_p,tau_t,tau_p,Tset,T,Pset,P)

%-------------------------------------------------------
%------------SIMULATION INITIALIZATION------------------
%-------------------------------------------------------
            
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

            
            
