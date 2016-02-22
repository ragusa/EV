function [flim,Wminus,Wplus] = compute_limited_flux_sums(...
    u_old,uFCT,dt,ML,AL,b,F,...
    sigma_min,sigma_max,source_min,source_max,theta,n_dof,...
    speed,inc,periodic_BC,limiting_option,DMP_option,...
    dirichlet_limiting_coefficient)            

% compute solution bounds
[Wplus,Wminus] = compute_DMP(u_old,uFCT,dt,ML,AL,b,theta,inc,periodic_BC);
if (DMP_option == 2) % max/min(DMP,CMP)
    [WplusCMP,WminusCMP] = compute_CMP(...
        u_old,sigma_min,sigma_max,source_min,source_max,speed*dt,...
        inc,periodic_BC);
    Wplus = max(Wplus,WplusCMP);
    Wminus = min(Wminus,WminusCMP);
end

% comput limited flux bounds
[Qplus,Qminus] = compute_Q(...
    u_old,uFCT,ML,Wplus,Wminus,AL,b,dt,theta);

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
