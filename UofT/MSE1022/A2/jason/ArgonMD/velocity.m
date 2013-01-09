function [v]=velocity(F,F_half,x,t_step,m);

a = F/m(1,1);
a_half = F_half/m(1,1);

for i=1:N
    for k=1:3
        if barostat==1   
        else
            v(i,k) = v(i,k) + t_step*(a(i,k)+a_half(i,k))/2 ;
        end
        
        if x(i,k)>=L(1,k)   %should it be <=???????
            x(i,k) = x(i,k) - L(1,k);
        end
        if x(i,k)<=0
            x(i,k) = x(i,k) + L(1,k);
        end
        
    end
end
            


% for i=1:N
%     for k=1:3
%         if barostat==1
%             x(i,k) = (1 + eps_p*t_step + t_step*t_step*(0.5/(tau_p*tau_p))*(P-Pset) + (0.5*t_step*t_step*eps_p*eps_p))*x(i,k) + p_half(i,k)*t_step/m(1,i) + 0.5*t_step*t_step*eps_p*p(i,k)/m(1,i);
%         else
%             x(i,k) = 
%         end
%         if x(i,k)>=L(1,k)   %should it be <=???????
%             x(i,k) = x(i,k) - L(1,k);
%         end
%         if x(i,k)<=0
%             x(i,k) = x(i,k) + L(1,k);
%         end
%     end
% end

