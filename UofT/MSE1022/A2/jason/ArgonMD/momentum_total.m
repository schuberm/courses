function p_sum=momentum_total(N,p)
p_sum(1,1:3)=0;
for k=1:3
    for i=1:N
        p_sum(1,k) = p_sum(1,k) + p(i,k);
    end
end




