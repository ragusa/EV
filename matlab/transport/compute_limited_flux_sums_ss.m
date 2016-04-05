function [flim,Wminus,Wplus] = compute_limited_flux_sums_ss(u,uL,F,AL_mod,b_mod,...
    sigma_min,sigma_max,source_min,source_max,mesh,phys,n_dof,fct_opts)

% unpack options
DMP_option = fct_opts.DMP_option;
limiting_option = fct_opts.limiting_option;
enforce_antidiffusion_bounds_signs = ...
    fct_opts.enforce_antidiffusion_bounds_signs;
dirichlet_limiting_coefficient = fct_opts.dirichlet_limiting_coefficient;

% compute max principle bounds
if (DMP_option == 1)
    [Wplus,Wminus] = compute_DMP_ss(u,AL_mod,b_mod,phys.inc);
elseif (DMP_option == 2)
    [Wplus,Wminus] = compute_DMP_ss(u,AL_mod,b_mod,phys.inc);
    [Wplus_analytic,Wminus_analytic] = compute_analytic_bounds_ss(...
        u,sigma_min,sigma_max,source_min,source_max,mesh.dx_min,phys.inc);
    Wplus  = max(Wplus, Wplus_analytic);
    Wminus = min(Wminus,Wminus_analytic);
elseif (DMP_option == 3)
    [Wplus,Wminus] = compute_analytic_bounds_ss(...
        u,sigma_min,sigma_max,source_min,source_max,mesh.dx_min,phys.inc,false);
elseif (DMP_option == 4)
    [Wplus,Wminus] = compute_analytic_bounds_ss(...
        u,sigma_min,sigma_max,source_min,source_max,mesh.dx_min,phys.inc,true);
else
    error('Invalid FCT solution bounds option');
end

% compute antidiffusion bounds
[Qplus,Qminus] = compute_Q_ss(u,Wplus,Wminus,AL_mod,b_mod);

% enforce signs of antidiffusion bounds if requested
if (enforce_antidiffusion_bounds_signs)
    Qplus  = max(Qplus, 0);
    Qminus = min(Qminus, 0);
    [Wplus,Wminus] = compute_W_from_Q_ss(u,Qplus,Qminus,AL_mod,b_mod);
end

% compute limiting coefficients
switch limiting_option
    case 0 % Full limiter
        flim = zeros(n_dof,1);
    case 1 % No limiter
        flim = sum(F,2);
    case 2 % Zalesak limiter
        flim = limiter_zalesak(F,Qplus,Qminus,phys.periodic_BC,...
            dirichlet_limiting_coefficient);
    case 3 % Josh limiter
        flim = limiter_josh(F,Qplus,Qminus,phys.periodic_BC,...
            dirichlet_limiting_coefficient);
    otherwise
        error('Invalid limiting option');
end

end
