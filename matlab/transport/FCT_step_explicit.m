function uFCT = FCT_step_explicit(u_old,uH,dt,ML,MC,AL,DH,DL,b,inc,speed,...
    sigma_min,sigma_max,q_min,q_max,DMP_option,limiting_option,...
    periodic_BC,modify_for_strong_DirichletBC,prelimit)

% size of system
n_dof = length(u_old);

% compute solution bounds
[Wplus,Wminus] = compute_DMP(u_old,u_old,dt,ML,AL,b,0,inc,periodic_BC);
if (DMP_option == 2) % max/min(DMP,CMP)
    [WplusCMP,WminusCMP] = compute_CMP(...
        u_old,sigma_min,sigma_max,q_min,q_max,speed*dt,inc,periodic_BC);
    Wplus = max(Wplus,WplusCMP);
    Wminus = min(Wminus,WminusCMP);
end

% theta = 0 for explicit
theta = 0;

% compute limited flux bounds
[Qplus,Qminus] = compute_Q(...
    u_old,u_old,ML,Wplus,Wminus,AL,b,dt,theta);

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
        flim = limiter_zalesak(F,Qplus,Qminus,periodic_BC);
    case 3 % Josh's limiter
        flim = limiter_josh(F,Qplus,Qminus,periodic_BC);
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
