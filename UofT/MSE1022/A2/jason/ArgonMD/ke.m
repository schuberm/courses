function kinetic=ke(N,p,m)
%calculate Kinetic Energy

kinetic=0;
for i=1:N
    kinetic = kinetic + 0.5*(1/m(1,i))*(p(i,1)*p(i,1) +p(i,2)*p(i,2) + p(i,3)*p(i,3));
end
    
