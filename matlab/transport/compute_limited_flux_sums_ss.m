function [flim,Wminus,Wplus] = compute_limited_flux_sums_ss(u,F,AL_mod,b_mod,...
    sigma_min,sigma_max,source_min,source_max,phys,n_dof,fct_opts)

% unpack options
DMP_option = fct_opts.DMP_option;
limiting_option = fct_opts.limiting_option;
dirichlet_limiting_coefficient = fct_opts.dirichlet_limiting_coefficient;

% compute max principle bounds
if DMP_option == 1
    [Wplus,Wminus] = compute_DMP_ss(u,AL_mod,b_mod,phys.inc);
elseif DMP_option == 2
    [Wplus,Wminus] = compute_DMP_ss(u,AL_mod,b_mod,phys.inc);
    [Wplus_analytic,Wminus_analytic] = compute_analytic_bounds_ss(...
        u,sigma_min,sigma_max,source_min,source_max,0,phys.inc);
    Wplus  = max(Wplus, Wplus_analytic);
    Wminus = min(Wminus,Wminus_analytic);
else
    error('Invalid FCT solution bounds option');
end

[Qplus,Qminus] = compute_Q_ss(u,Wplus,Wminus,AL_mod,b_mod);

% compute limiting coefficients
switch limiting_option
    case 0 % Full limiter
        flim = zeros(n_dof,1);
    case 1 % No limiter
        flim = sum(F,2);
    case 2 % Zalesak limiter
        flim = limiter_zalesak(F,Qplus,Qminus,phys.periodic_BC,
            dirichlet_limiting_coefficient);
    case 3 % Josh limiter
        flim = limiter_josh(F,Qplus,Qminus,phys.periodic_BC,
            dirichlet_limiting_coefficient);
    otherwise
        error('Invalid limiting option');
end

end
