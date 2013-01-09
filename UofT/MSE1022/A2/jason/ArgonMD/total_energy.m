function [KE,PE,TE]=total_energy(N,x,L,cutoff,f_cutoff,a,pe_cutoff,v,m,a2,epsilon_Ar,sigma_Ar)

KE = kineticenergy(N,v,m);
PE = 0;
PE=Pot_Energy(N,L,pe_cutoff,a,a2,f_cutoff,cutoff,x,epsilon_Ar,sigma_Ar)
TE = KE + PE;


% C Code
% 			// output 
% 			if (t%stats == 0) {			
% 				// calculate PE, KE and TE
% 				KE = ke(p,m);
% 				T = (2.*ke(p,m))/(3.*(N-1.));
% 				PE = 0.;
% 				for (i=0; i<N; i++){
% 					for (j=i+1; j<N; j++){
% 						for (k=0; k<3; k++) {
% 							rij[k] = x[i][k] - x[j][k];
% 							if (rij[k] > (L[k]/2.)) rij[k] = rij[k] - L[k];    // nearest image
% 							if (rij[k] < (-L[k]/2.)) rij[k] = rij[k] + L[k];
% 						}
% 						r2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
% 						if (r2 < a2) {
% 							if (cutoff == true) PE = PE + pe(r2) - pe_cutoff;
% 							else {
% 								PE = PE + pe(r2) - pe_cutoff + (f_cutoff*a*(sqrt(r2) - a));	
% 							}
% 						}
% 					}	
% 				}
% 				TE = KE + PE;


