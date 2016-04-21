function [flim,Wminus,Wplus] = compute_limited_flux_sums(...
    u_old,uFCT,dt,ML,AL,b,F,...
    sigma_min,sigma_max,source_min,source_max,mesh,theta,n_dof,...
    phys,fct_opts)

% unpack options
speed = phys.speed;
inc = phys.inc;
periodic_BC = phys.periodic_BC;
limiting_option = fct_opts.limiting_option;
DMP_option = fct_opts.DMP_option;
enforce_antidiffusion_bounds_signs = ...
    fct_opts.enforce_antidiffusion_bounds_signs;
dirichlet_limiting_coefficient = fct_opts.dirichlet_limiting_coefficient;

% compute solution bounds
if (DMP_option == 1)
    [Wplus,Wminus] = compute_DMP(u_old,uFCT,dt,ML,AL,b,theta,inc,periodic_BC);
elseif (DMP_option == 2)
    [Wplus,Wminus] = compute_DMP(u_old,uFCT,dt,ML,AL,b,theta,inc,periodic_BC);
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

% comput limited flux bounds
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

% compute limited fluxes
switch limiting_option
    case 0 % no correction
        flim = zeros(n_dof,1);
    case 1 % full correction (no limiting)
        flim = sum(F,2);
    case 2 % Zalesak's limiter
        flim = limiter_zalesak(F,Qplus,Qminus,periodic_BC,...
            dirichlet_limiting_coefficient);
    case 3 % Josh's limiter
        flim = limiter_josh(F,Qplus,Qminus,periodic_BC,...
            dirichlet_limiting_coefficient);
    otherwise
        error('Invalid limiting option');
end

end
