function p_half=half_momentum(N,p,F,t_step,eta_t,eps_p,p_half)
%p_half(1:N,1:3)=zeros(N,3); 

for i=1:N
    for k=1:3
            p_half(i,k) = p(i,k) + (F(i,k) - (eta_t+eps_p)*p(i,k))*t_step*0.5;
        end
    end
end

