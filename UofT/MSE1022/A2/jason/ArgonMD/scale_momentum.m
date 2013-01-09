function p=scale_momentum(N,p)

    % scale initial momenta to zero
            p_sum(1,1:3)=0; 
            for k=1:3
            p_sum(1,k) = p_sum(1,k)+sum(p(:,k)); % sum momentum
            end
            p_sub(1,1:3) = p_sum(:,1:3)/N;	% find amount to subtract from each atom
            for k=1:3
            p(:,k) = p(:,k)-p_sub(1,k); % subtract momentum
            end

