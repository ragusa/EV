function [flim,Wminus,Wplus] = compute_limited_flux_sums_ss(u,F,AL_mod,b_mod,...
    sigma_min,sigma_max,source_min,source_max,mesh,phys,n_dof,fct_opts,opts,...
    source_over_sigma_min, source_over_sigma_max)

% unpack options
DMP_option = fct_opts.DMP_option;
enforce_antidiffusion_bounds_signs = ...
    fct_opts.enforce_antidiffusion_bounds_signs;
dirichlet_limiting_coefficient = fct_opts.dirichlet_limiting_coefficient;

% compute max principle bounds
if (DMP_option == 1)
    [Wplus,Wminus] = compute_DMP_ss(u,AL_mod,b_mod,phys.inc,opts);
elseif (DMP_option == 2)
    [Wplus,Wminus] = compute_DMP_ss(u,AL_mod,b_mod,phys.inc,opts);
    [Wplus_analytic,Wminus_analytic] = compute_analytic_bounds_ss(...
        u,sigma_min,sigma_max,source_min,source_max,mesh.dx_min,phys.inc,false);
    Wplus  = max(Wplus, Wplus_analytic);
    Wminus = min(Wminus,Wminus_analytic);
elseif (DMP_option == 3)
    [Wplus,Wminus] = compute_analytic_bounds_ss(...
        u,sigma_min,sigma_max,source_min,source_max,mesh.dx_min,phys.inc,false);
elseif (DMP_option == 4)
    [Wplus,Wminus] = compute_analytic_bounds_ss(...
        u,sigma_min,sigma_max,source_min,source_max,mesh.dx_min,phys.inc,true);
elseif (DMP_option == 5)
    [Wplus,Wminus] = compute_analytic_bounds_ss_alternate(...
        u,sigma_min,sigma_max,source_over_sigma_min,source_over_sigma_max,mesh.dx_min,phys.inc,true);
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
