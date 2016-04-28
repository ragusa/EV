function [uFCT,Wminus,Wplus] = FCT_step_explicit(u_old,uH,dt,ML,MC,AL,DH,DL,b,inc,speed,...
    sigma_min,sigma_max,q_min,q_max,mesh,fct_opts,opts,periodic_BC)

% unpack options
DMP_option = fct_opts.DMP_option;
limiting_option = fct_opts.limiting_option;
prelimit = fct_opts.prelimit;
dirichlet_limiting_coefficient = fct_opts.dirichlet_limiting_coefficient;
modify_for_strong_DirichletBC = opts.modify_for_strong_DirichletBC;

% size of system
n_dof = length(u_old);

% compute solution bounds
if (DMP_option == 1)
    [Wplus,Wminus] = compute_DMP(u_old,u_old,dt,ML,AL,b,0,inc,periodic_BC);
elseif (DMP_option == 2)
    [Wplus,Wminus] = compute_DMP(u_old,u_old,dt,ML,AL,b,0,inc,periodic_BC);
    [WplusCMP,WminusCMP] = compute_analytic_bounds(...
        u_old,sigma_min,sigma_max,q_min,q_max,speed*dt,mesh,inc,periodic_BC);
    Wplus = max(Wplus,WplusCMP);
    Wminus = min(Wminus,WminusCMP);
else
    error('Invalid FCT solution bounds option');
end

% theta = 0 for explicit
theta = 0;

% compute limited flux bounds
[Qplus,Qminus] = compute_Q(...
    u_old,u_old,ML,Wplus,Wminus,AL,b,dt,theta);

% check signs of antidiffusion bounds Q
%check_antidiffusion_bounds_signs(Qplus,Qminus);

% compute flux corrections
F = flux_correction_matrix(u_old,uH,dt,DH,DL,MC,theta);

% prelimit flux corrections if user specified
if (prelimit)
    % compute low-order solution
    uL = compute_low_order_solution_theta(u_old,AL,...
        ML,b,dt,0,inc,modify_for_strong_DirichletBC);

    % prelimit fluxes
    F = prelimit_fluxes(F,uL);
end
        
% compute limited fluxes
switch limiting_option
    case 0 % no correction
        flim = zeros(n_dof,1);
    case 1 % full correction (no limiting)
        flim = sum(F,2);
    case 2 % Zalesak's limiter
        flim = limiter_zalesak(F,Qplus,Qminus,opts,dirichlet_limiting_coefficient);
    case 3 % Josh's limiter
        flim = limiter_josh(F,Qplus,Qminus,opts,dirichlet_limiting_coefficient);
    otherwise
        error('Invalid limiting option');
end

% create system matrix and rhs
ML_mod = ML;
rhs = (ML - dt*AL)*u_old + dt*(b + flim);
if modify_for_strong_DirichletBC
    ML_mod(1,:)=0; ML_mod(1,1)=1;
    rhs(1) = inc;
end

% compute FCT solution
uFCT = ML_mod \ rhs;

end
