function [x,x_p]=position(N,t_step,x,x_p,p_half,p,m,L,tau_p,eps_p,P,Pset)

term1 = eps_p*t_step;
term2 = (t_step*t_step/(2*(tau_p*tau_p)))*(P-Pset);
term3 = (t_step*eps_p)*(t_step*eps_p)*0.5;

for i=1:N
    for k=1:3
        x(i,k) = (1 + term1 + term2 + term3)*x(i,k) + p_half(i,k)*t_step/m(1,i) + (t_step*t_step)*eps_p*p(i,k)*0.5/m(1,i);
        x_p(i,k) = (1 + term1 + term2 + term3)*x_p(i,k) + (t_step*t_step)*eps_p*p(i,k)*0.5/m(1,i);
        if x(i,k)>L(1,k)   %should it be <=???????
            x(i,k) = x(i,k) - L(1,k);
        end
        if x(i,k)<0
            x(i,k) = x(i,k) + L(1,k);
        end
        
    end
end
            

