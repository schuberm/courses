function part_PHI=partPE(N,part,x)
%calculate Potential Energy of one particle

    part_PHI = 0;
    
    for i=1:N
        if i~=part
            for k=1:3
                rij(1,k) = x(part,k) - x(i,k);
                rij(1,k) = periodic(rij(1,k),L(1,k));   %periodic boundaries
            end
            r2 = rij(1,1)*rij(1,1) +rij(1,2)*rij(1,2) + rij(1,3),*rij(1,3);
            if r2 < a2
                if cutoff==1
                    part_PHI = part_PHI +0.5*(pe(r2) - pe_cutoff);
                else
                    part_PHI = part_PHI + 0.5*(pe(r2) - pe_cutoff +(f_cutoff*a*(sqrt(r2)-a)));
                end
            end
        end
    end
    
