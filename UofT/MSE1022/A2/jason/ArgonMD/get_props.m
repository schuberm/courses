function [props,pIndex]=get_props(x,L,t,t_stats,N,cutoff,f_cutoff,a,pe_cutoff,v,m,a2,epsilon_Ar,sigma_Ar,pIndex,props,kb)

            if rem(t,t_stats)==0 
            plot_pos2D(x,L);
            [KE,PE,TE]=total_energy(N,x,L,cutoff,f_cutoff,a,pe_cutoff,v,m,a2,epsilon_Ar,sigma_Ar)
            T=temperature(KE,N,kb)
            props(pIndex,1)=t;
            props(pIndex,2)=T;
            props(pIndex,3)=KE;
            props(pIndex,4)=PE;
            props(pIndex,5)=TE;
            pIndex=pIndex+1;
            
            sumsq=v.*v;
            sumv2(:,1)=sumsq(:,1); sumv2(N+1:2*N,1)=sumsq(:,2); sumv2(2*N+1:3*N,1)=sumsq(:,3);
            [Prob, bin]=hist(sumsq(:,1),50); Prob=Prob/sum(Prob);
            loglog(bin,Prob,bin,bin.^(-1)/(bin(1,1)^(-1)))
            end

