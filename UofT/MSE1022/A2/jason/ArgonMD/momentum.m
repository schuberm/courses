function p=momentum(N,p_half,F,t_step,eta_t,eps_p)
p(1:N,1:3)=zeros(N,3); 
bottom = 1 + (eta_t + eps_p)*t_step*0.5;

for i=1:N
    for k=1:3
            p(i,k) = (p_half(i,k) + F(i,k)*0.5*t_step)/bottom;
        end
    end
end



