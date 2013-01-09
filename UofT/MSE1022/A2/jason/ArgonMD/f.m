function f=force_eval(r2,epsilon_Ar,sigma_Ar)

	r4 = r2*r2;
	r8 = r4*r4;
	inverse_r14 = 1./(r8*r4*r2);

	f = 48.*inverse_r14 - 24.*r4*r2*inverse_r14;

Fold=F;
%clear old Forces
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
                f = force(r2,epsilon_Ar,sigma_Ar) - f_cutoff*a*(1/sqrt(r2));
            end
            for k=1:3    
                F(i,k) = F(i,k) + rij(1,k)*f;
                F(j,k) = F(j,k) - rij(1,k)*f;
                v(i,k) = v(i,k) + 0.5*(Fold(i,k)+F(i,k))*t_step;
                v(j,k) = v(j,k) + 0.5*(Fold(j,k)+F(j,k))*t_step; 
            end
            P_viral = P_viral +r2*f;
        end
    end
end
P = (P_viral/(3.*V)) + ((N*T)/V);
%FIXXXXX!!!! http://www.sklogwiki.org/SklogWiki/index.php/Virial_pressure
% C Code
% 			// find forces/pressure at t+delt
% 			for (i=0; i<N; i++){					// clear old values
% 				for (k=0; k<3; k++) F[i][k] = 0.;
% 			}
% 			P_viral = 0.;		
% 			for (i=0; i<N; i++){
% 				for (j=i+1; j<N; j++){
% 					for (k=0; k<3; k++) {
% 						rij[k] = x[i][k] - x[j][k];
% 						if (rij[k] > (0.5*L[k])) rij[k] = rij[k] - L[k];    // nearest image
% 						if (rij[k] < (-0.5*L[k])) rij[k] = rij[k] + L[k];
% 					}
% 					r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
% 					if (r2 < a2) {
% 						if (cutoff == true) f = force(r2);
% 						else f = force(r2) - f_cutoff*a*(1./(sqrt(r2)));
% 						for (k=0; k<3; k++){
% 							F[i][k] = F[i][k] + rij[k]*f;
% 							F[j][k] = F[j][k] - rij[k]*f;
% 						}
% 						P_viral = P_viral + r2*f;
% 					}
% 				}
% 			}


