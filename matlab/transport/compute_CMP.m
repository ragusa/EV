function [Wplus,Wminus] = compute_CMP(u,sigma_min,sigma_max,q_min,q_max,s,inc,periodic_BC)

% size of system
n = length(u);

Wplus  = zeros(n,1);
Wminus = zeros(n,1);

% compute bounds
for i = 1:n
    % compute index range of support of i
    if ~periodic_BC
        iL = max(i-1,1);
        iR = min(i+1,n);
    else
        if i == 1
            iL = n;
            iR = 2;
        elseif i == n
            iL = n-1;
            iR = 1;
        else
            iL = i-1;
            iR = i+1;
        end
    end
    % compute max and min old solution in the support of each dof
    support = [iL i iR];
    u_max = max(u(support));
    u_min = min(u(support));
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

if ~periodic_BC
    Wplus(1) = inc;
    Wminus(1) = inc;
end

end