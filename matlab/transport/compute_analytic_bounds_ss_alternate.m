function [Wplus,Wminus] = compute_analytic_bounds_ss_alternate(...
    u,sigma_min,sigma_max,q_over_sigma_min,q_over_sigma_max,s,inc,upwind_only)

% size of system
n = length(u);

Wplus  = zeros(n,1);
Wminus = zeros(n,1);

% compute bounds
Wplus(1)  = inc;
Wminus(1) = inc;
for i = 2:n
    % get set of neighbor indices
    if (upwind_only)
        % upwind range, assuming flow is in +x direction
        i1 = max(i-1,1);
        i2 = min(i,n);
    else
        i1 = max(i-1,1);
        i2 = min(i+1,n);
    end
    % find min and max in support of i
    u_max = max(u(i1:i2));
    u_min = min(u(i1:i2));
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
