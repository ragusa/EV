function [Wplus,Wminus] = compute_analytic_bounds_alternate(...
   u,sigma_min,sigma_max,q_over_sigma_min,q_over_sigma_max,...
   s,mesh,inc,periodic_BC)

% size of system
n = length(u);

Wplus  = zeros(n,1);
Wminus = zeros(n,1);

% compute range of min/max operations
range = ceil(s/mesh.dx(1));

% compute bounds
Wplus(1)  = inc;
Wminus(1) = inc;
for i = 2:n
    % compute index range of support of i
    if (periodic_BC)
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
    else
        iL = max(i-range,1);
        iR = min(i+range,n);
    end
    % compute max and min old solution in the support of each dof
    u_max = max(u(iL:iR));
    u_min = min(u(iL:iR));
    % compute source bounds
    if (sigma_min(i) == 0) % assume no vacuum
        error('No vacuum allowed for these solution bounds');
    end
    src_max = q_over_sigma_max(i)*(1 - exp(-sigma_min(i)*s));
    src_min = q_over_sigma_min(i)*(1 - exp(-sigma_max(i)*s));

    % compute upper and lower bounds
    Wplus(i)  = u_max*exp(-sigma_min(i)*s) + src_max;
    Wminus(i) = u_min*exp(-sigma_max(i)*s) + src_min;
end

end
