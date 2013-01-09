function out = fcclattice(l,m,n)
%
% This function creates a list of coordinates of fcc lattice points.
% l, m, and n are the number of unit cells in the x, y, and z-direction, 
% The lattice constant is 1.
%

N = l*m*n;
[x,y,z] = meshgrid(1:l,1:m,1:n);
p(1:N,1) = x(:);  p(N+1:2*N,1) = x(:)+0.5;  p(2*N+1:3*N,1) = x(:)+0.5;  p(3*N+1:4*N,1) = x(:);
p(1:N,2) = y(:);  p(N+1:2*N,2) = y(:)+0.5;  p(2*N+1:3*N,2) = y(:);  p(3*N+1:4*N,2) = y(:)+0.5;
p(1:N,3) = z(:);  p(N+1:2*N,3) = z(:);  p(2*N+1:3*N,3) = z(:)+0.5;  p(3*N+1:4*N,3) = z(:)+0.5;
out = p;