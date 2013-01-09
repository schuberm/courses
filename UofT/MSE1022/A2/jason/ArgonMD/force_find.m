function [F,P,P_viral]=force_find(N,L,a2,f_cutoff,a,x,cutoff,V,T)
%clear old Forces
F(1:N,1:3)=zeros(N,3);
% rij(1,1:3) = 

P_viral = 0;
for i=1:N
%     rij(1:N-i,1:3) = 
        for k=1:3
            rij(:,k) = x(i,k) - x(i+1:N,k);
            rij(:,k) = periodic(rij(:,k),L(:,k));   %periodic boundaries
        end
        r2 = (x(i,1)-x(i+1:N,1)).^2 + (x(i,2) - x(i+1:N,2)).^2 + (x(i,3) - x(i+1:N,3)).^2;         
        [I]=find(r2<a2);
        %[J]=find(I~=i);
            if cutoff == 1
                f = force_eval_find(r2(I));
            else
                f = force_eval_find(r2(I));
            end
            for k=1:3    
                F(i,k) = F(i,k) + sum(rij(:,k).*f);
                F(i+1:N,k) = F(i+1:N,k) - rij(:,k).*f;
            end
            P_viral = P_viral +sum(r2.*f);
            clear rij
end
         
%P = (P_viral/(3*V)) + ((N*T)/V);
%Set P and Pset to 1 for the purpose of NanoMD
P=1;



