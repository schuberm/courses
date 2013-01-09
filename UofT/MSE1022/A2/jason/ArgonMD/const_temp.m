function v=const_temp(N,v,m,Tset,kb)

%I use this for single species.  Multiple species would need v=p./m
mass=m(1,1);
tfac = 3*N*kb*Tset/mass;
sumvsq=sum(sum(v.*v));
fac = sqrt(tfac/sumvsq);
v=v*fac;

% KE = ke(N,p,m);
% KE_set = 1.5*(N-1)*Tset;
% alpha = (KE_set)/KE;
% alpha = sqrt(alpha);
% p = p*alpha;

% C Code
% he has this commented out
% 			// scale momentum to desired T (constant kinetic energy ensemble)
% 			// specify the time range where this is to be applied
% 	/*		if (t >= 0 && t < 10000) {
% 				KE = ke(p,m);
% 				KE_set = 1.5*(N-1)*Tset;	
% 				alpha = (KE_set)/KE;
% 				alpha = sqrt(alpha);
% 				for (i=0; i<N; i++){
% 					for (k=0; k<3; k++) p[i][k] = p[i][k]*alpha;
% 				}
% 			} */

