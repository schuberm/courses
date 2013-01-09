function kinetic=kineticenergy(N,v,m)
%calculate Kinetic Energy

kinetic=0;
for i=1:N
    kinetic = kinetic + 0.5*m(1,i)*(v(i,1)*v(i,1) +v(i,2)*v(i,2) + v(i,3)*v(i,3));
end
    
