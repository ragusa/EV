function [Wplus,Wminus] = compute_CMP_ss(u,sigma_min,sigma_max,q_min,q_max,s,inc)

% size of system
n = length(u);

Wplus  = zeros(n,1);
Wminus = zeros(n,1);

% compute bounds
Wplus(1)  = inc;
Wminus(1) = inc;
for i = 2:n
    % get set of neighbor indices
    i1 = max(i-1,1);
    i2 = min(i+1,n);
    % find min and max in support of i
    u_max = max(u(i1:i2));
    u_min = min(u(i1:i2));
    % compute source bounds
    if (sigma_min(i) == 0)
        src_max = q_max(i)*s*exp(-sigma_min(i)*s);
    else
        src_max = q_max(i)/sigma_min(i)*(1 - exp(-sigma_min(i)*s));
    end
    if (sigma_max(i) == 0)
        src_min = q_min(i)*s*exp(-sigma_max(i)*s);
    else
        src_min = q_min(i)/sigma_max(i)*(1 - exp(-sigma_max(i)*s));
    end
    % compute upper and lower bounds
    Wplus(i)  = u_max*exp(-sigma_min(i)*s) + src_max;
    Wminus(i) = u_min*exp(-sigma_max(i)*s) + src_min;
end

end
