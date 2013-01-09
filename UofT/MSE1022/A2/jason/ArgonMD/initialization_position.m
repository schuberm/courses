function [x,V,L]=initialization_position(N,L,m,mass_Ar,a_0,ncell)

%-------------------------------------------------------
%------------SIMULATION INITIALIZATION------------------
%-------------------------------------------------------

%Produce Intial lattice
out = fcclattice(ncell,ncell,ncell); x = a_0*(out-1);
%Set initial species mass
m(1,1:N) = mass_Ar; 
%Randomize Initial Coordiantes away from Equillirbium by 5% a_0
% latper = 0.05;
% x = x +a_0*latper*(rand(N,3)-0.5);
% %Adjust the positions so that all random coords are inside box
% x = x + abs(min(min(x)))+a_0*latper;
% x_0=x;
% simulation cell size
	for k=1:3 
        L(1,k) = max(max(x))+a_0/2;	
        % desired system size. May be specified here, determined from input coords.
        %So, I add an extra a_0 so that using the periodic boundary
        %conditions adds an extra lattice constant in between mirrored
        %particles.
    end
	V = L(1,1)*L(1,2)*L(1,3);    %system volume
    plot_pos(x,L);


            
            
