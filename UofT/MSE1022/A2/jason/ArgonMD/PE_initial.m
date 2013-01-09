function PE=PE_initial(N,PE,rij,L,pe_cutoff,a,a2,f_cutoff,cutoff,x)


for i=1:N
    for j=(i+1):N
        for k=1:3
            rij(1,k) = x(i,k) - x(j,k);
            rij(1,k) = periodic(rij(1,k),L(1,k));   %periodic boundaries
        end
        r2 = rij(1,1)^2 + rij(1,2)^2 + rij(1,3)^2;
        
        if (r2<a2)
            if cutoff==1
                PE = PE+pe(r2) - pe_cutoff;
            else
                PE = PE + pe_cutoff + (f_cutoff*a*(sqrt(r2) - a));
            end
        end
    end
end



% C Code
% 		for (i=0; i<N; i++){
% 			for (j=i+1; j<N; j++){
% 				for (k=0; k<3; k++) {
% 					rij[k] = x[i][k] - x[j][k];
% 					if (rij[k] > (L[k]/2.)) rij[k] = rij[k] - L[k];    % nearest image
% 					if (rij[k] < (-L[k]/2.)) rij[k] = rij[k] + L[k];
% 				}
% 				r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
% 				if (r2 < a2) {
% 					if (cutoff == true) PE = PE + pe(r2) - pe_cutoff;
% 					else {
% 						PE = PE + pe(r2) - pe_cutoff + (f_cutoff*a*(sqrt(r2) - a));	
% 					}
% 				}	
% 			}	
% 		}
