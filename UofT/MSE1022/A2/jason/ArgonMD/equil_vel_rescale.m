function [x,v,F,P,P_viral,T,KE,PE]=equil_vel_rescale(t_total,N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x,epsilon_Ar,sigma_Ar,m,kb,V,t_step,mass_Ar)

    %-------------------------------------------------------
    %------------EQUILIBRIATION-----------------------------
    %-------------------------------------------------------
    %Run Equilibration at constant T to equilibriate.
            for t=1:t_total
                [x]=position(N,barostat,t_step,x,F,v,m,L);
                % find forces/pressure/velocity at t+delt
                [F,P,P_viral,v]=force(N,F,L,a2,f_cutoff,a,x,cutoff,epsilon_Ar,sigma_Ar,v,t_step,V,T,m);
                %Use simple velocity rescaling
                v=vel_rescale(N,v,m,Tset,kb);
                [props,pIndex]=get_props(x,L,t,t_stats,N,cutoff,f_cutoff,a,pe_cutoff,v,m,a2,epsilon_Ar,sigma_Ar,pIndex,props,kb);
            end 
