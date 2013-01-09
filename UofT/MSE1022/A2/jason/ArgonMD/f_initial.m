function [F,P_viral,P]=f_initial(cutoff,N,x,L,a,a2,V,T,epsilon_Ar,sigma_Ar,f_cutoff)

F(1:N,1:3)=0;   %clear old values
P_viral = 0;
for i=1:N
    for j=(i+1):N
        for k=1:3
            rij(1,k) = x(i,k) - x(j,k);
            rij(1,k) = periodic(rij(1,k),L(1,k));   %periodic boundaries
        end
            r2 = rij(1,1)^2 + rij(1,2)^2 + rij(1,3)^2;
            if r2<a2
                if cutoff == 1
                    f = force(r2,epsilon_Ar,sigma_Ar);
                else
                    f=force(r2,epsilon_Ar,sigma_Ar) - f_cutoff*a*(1/sqrt(r2));
                end
                for k=1:3
                F(i,k) = F(i,k) + rij(1,k)*f;
                F(j,k) = F(j,k) - rij(1,k)*f;
                end
               P_viral = P_viral +r2*f; 
            end
    end
end
P = (P_viral/(3.*V)) + ((N*T)/V);
%FIXXXXX!!!! http://www.sklogwiki.org/SklogWiki/index.php/Virial_pressure



%C code
% double f;
% 		for (i=0; i<N; i++){					// clear old values
% 			for (k=0; k<3; k++) F[i][k] = 0.;
% 		}
% 		for (i=0; i<N; i++){
% 			for (j=i+1; j<N; j++){
% 				for (k=0; k<3; k++) {
% 					rij[k] = x[i][k] - x[j][k];
% 					if (rij[k] > (L[k]/2.)) rij[k] = rij[k] - L[k];    // nearest image
% 					if (rij[k] < (-L[k]/2.)) rij[k] = rij[k] + L[k];
% 				}
% 				r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
% 				if (r2 < a2) {
% 					if (cutoff == true) f = force(r2);
% 					else f = force(r2) - f_cutoff*a*(1./(sqrt(r2)));
% 					for (k=0; k<3; k++){
% 						F[i][k] = F[i][k] + rij[k]*f;
% 						F[j][k] = F[j][k] - rij[k]*f;
% 					}
% 					P_viral = P_viral + r2*f;
% 				}
% 			}
% 		}
% P = (P_viral/(3.*V)) + ((N*T)/V);