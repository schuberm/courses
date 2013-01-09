function [x,p,F,P,P_viral,T,KE,PE]=initialization_vel(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,m,x,Tset,V)

%-------------------------------------------------------
%------------SIMULATION INITIALIZATION------------------
%-------------------------------------------------------

            PE=Pot_Energy(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x);
            %Set Initial velocities Randomly
            p = initial_vel(N);
            p=vel_rescale(N,p,m,Tset);
%             for i=1:N
%             for k=1:3
%             p(i,k) = v(i,k).*m(1,i);
%             end
%             end
            KE = ke(N,p,m);
            T=temperature(N,p,m);
            % find initial forces/pressure
            [F,P,P_viral]=force(N,L,a2,f_cutoff,a,x,cutoff,V,T);
            % find initial barostat and thermostat parameters
            %baro_initial();


            
            
