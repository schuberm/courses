function xrms = mean_square(x_p,x_0,N)

for i=1:N
        xrms(i,1) = sqrt( (x_p(i,1)-x_0(i,1))^2 + (x_p(i,2)-x_0(i,2))^2 + (x_p(i,3)-x_0(i,3))^2);
end
xrms = mean(xrms);



