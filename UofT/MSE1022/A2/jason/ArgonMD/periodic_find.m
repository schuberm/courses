function r=periodic(r,L)
%determine periodic boundary conditions

[I]=find(r>(L/2));
r(I,1) = r(I,1) - L; r(I,2) = r(I,2) - L; r(I,3) = r(I,3) - L;
[J]=find(r<(-L/2));
r(I,1) = r(I,1) + L; r(I,2) = r(I,2) + L; r(I,3) = r(I,3) + L;


