function [flim,Wminus,Wplus] = compute_limited_flux_sums(...
    u_old,uFCT,uL,dt,ML,AL,b,F,...
    sigma_min,sigma_max,source_min,source_max,mesh,theta,n_dof,...
    phys,fct_opts,opts)

% unpack options
speed = phys.speed;
inc = phys.inc;
periodic_BC = phys.periodic_BC;
DMP_option = fct_opts.DMP_option;
enforce_antidiffusion_bounds_signs = ...
    fct_opts.enforce_antidiffusion_bounds_signs;
dirichlet_limiting_coefficient = fct_opts.dirichlet_limiting_coefficient;

% compute solution bounds
if (DMP_option == 1)
    [Wplus,Wminus] = compute_DMP(u_old,uL,dt,ML,AL,b,theta,inc,periodic_BC);
elseif (DMP_option == 2)
    [Wplus,Wminus] = compute_DMP(u_old,uL,dt,ML,AL,b,theta,inc,periodic_BC);
    [WplusCMP,WminusCMP] = compute_analytic_bounds(...
        u_old,sigma_min,sigma_max,source_min,source_max,speed*dt,mesh,...
        inc,periodic_BC);
    Wplus = max(Wplus,WplusCMP);
    Wminus = min(Wminus,WminusCMP);
elseif (DMP_option == 3)
    [Wplus,Wminus] = compute_analytic_bounds(...
        u_old,sigma_min,sigma_max,source_min,source_max,speed*dt,mesh,...
        inc,periodic_BC);
else
    error('Invalid FCT solution bounds option');
end

% check if solution bounds are already satisfied
%bounds_satisfied = check_solution_bounds(uFCT,Wplus,Wminus);

% compute limited flux bounds
[Qplus,Qminus] = compute_Q(...
    u_old,uFCT,ML,Wplus,Wminus,AL,b,dt,theta);

% enforce signs of antidiffusion bounds if requested
if (enforce_antidiffusion_bounds_signs)
    Qplus  = max(Qplus, 0);
    Qminus = min(Qminus, 0);
    [Wplus,Wminus] = compute_W_from_Q(u_old,uFCT,ML,Qplus,Qminus,AL,b,dt,theta);
end
if ~periodic_BC
    Wplus(1) = inc;
    Wminus(1) = inc;
end

% check signs of antidiffusion bounds Q
check_antidiffusion_bounds_signs(Qplus,Qminus);

% compute limited antidiffusion sums
if (fct_opts.use_multipass_limiting)
    flim = multipass_limiter(F,Qplus,Qminus,opts,fct_opts);
else
    [flim,~] = fct_opts.limiter(F,Qplus,Qminus,zeros(n_dof,1),...
        opts,dirichlet_limiting_coefficient);
end

end
