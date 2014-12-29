function [Wplus,Wminus] = compute_max_principle_bounds(u,AL,b,inc)

n = length(u);

Wplus = zeros(n,1);
Wminus = zeros(n,1);
Wplus(1)  = inc;
Wminus(1) = inc;
for i = 2:n
    % compute index range of support of i
    i1 = max(i-1,1);
    i2 = min(i+1,n);
    % compute max and min old solution in the support of each dof
    u_max = max(u(i1:i2));
    u_min = min(u(i1:i2));
    % compute bounds
    Wplus(i)  = (1-sum(AL(i,:))/AL(i,i))*u_max + b(i)/AL(i,i);
    Wminus(i) = (1-sum(AL(i,:))/AL(i,i))*u_min + b(i)/AL(i,i);
end

end