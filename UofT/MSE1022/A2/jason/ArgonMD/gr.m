function G=gr(N,x,L)
%calculate the Radial Distribution function G(r)
num=5000;
delr = L(1,1)/num;
rij(1,1:3)=0;
G(:,1)=(1:num)*delr;
G(:,2)=zeros(1,num);
rho = N/(L(1,1)*L(1,2)*L(1,3));
    
    for i=1:N
        for j=1:N
            if i~=j %don't count the particle itself
                for k=1:3
                    rij(1,k) = x(i,k) - x(j,k);
                    rij(1,k) = periodic(rij(1,k),L(1,k));   %periodic boundaries
                        r2 = rij(1,1)*rij(1,1) + rij(1,2)*rij(1,2) + rij(1,3)*rij(1,3);
                        r = sqrt(r2);
                        bin = floor(r/delr);
                        G(bin,2) = G(bin,2)+1;
                end
            end
        end
    end
    
    for i=1:num
        G(i,2) = G(i,2)/N;
        G(i,2) = G(i,2)*(1/(rho*4*pi*(delr*delr*delr*(i*i +i +0.25))));
    end
    
    %UNFINISHED: need to complete the output file
    
    %dlmwrite(
    
    
