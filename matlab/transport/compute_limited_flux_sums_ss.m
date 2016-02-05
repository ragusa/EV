function [flim,Wminus,Wplus] = compute_limited_flux_sums_ss(u,F,AL_mod,b_mod,...
    sigma_min,sigma_max,source_min,source_max,phys,n_dof,DMP_option,...
    limiting_option)

% compute max principle bounds
[Wplus,Wminus] = compute_DMP_ss(u,AL_mod,b_mod,phys.inc);
if DMP_option == 2
    [WplusCMP,WminusCMP] = compute_CMP_ss(u,sigma_min,...
        sigma_max,source_min,source_max,0,phys.inc);
    Wplus = max(Wplus,WplusCMP);
    Wminus = min(Wminus,WminusCMP);
end

[Qplus,Qminus] = compute_Q_ss(u,Wplus,Wminus,AL_mod,b_mod);

% compute limiting coefficients
switch limiting_option
    case 0 % Full limiter
        flim = zeros(n_dof,1);
    case 1 % No limiter
        flim = sum(F,2);
    case 2 % Zalesak limiter
        flim = limiter_zalesak(F,Qplus,Qminus,phys.periodic_BC);
    case 3 % Josh limiter
        flim = limiter_josh(F,Qplus,Qminus,phys.periodic_BC);
    otherwise
        error('Invalid limiting option');
end

end
