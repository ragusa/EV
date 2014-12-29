function uFCT = FCT_step_explicit(u_old,uH,dt,ML,MC,AL,DH,DL,b,inc,speed,...
    sigma_min,sigma_max,q_min,q_max,DMP_option,limiting_option,periodic_BC)

% size of system
n_dof = length(u_old);

% compute solution bounds
switch DMP_option
    case 1 % DMP
        [Wplus,Wminus] = compute_DMP(...
            u_old,u_old,dt,ML,AL,b,0,inc,periodic_BC);
    case 2 % CMP
        [Wplus,Wminus] = compute_CMP(...
            u_old,sigma_min,sigma_max,q_min,q_max,speed*dt,inc,periodic_BC);
    otherwise
        error('Invalid DMP option');
end
% compute limited flux bounds
[Qplus,Qminus] = compute_Q(...
    u_old,u_old,ML,Wplus,Wminus,AL,b,dt,0);

% compute flux corrections
F = flux_correction_matrix(u_old,uH,dt,DH,DL,MC,0);

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

% compute FCT solution
rhs = (ML - dt*AL)*u_old + dt*b + flim;
ML_mod = ML;
% modify matrix and rhs if Dirichlet BC are used
if ~periodic_BC
    rhs(1) = inc;
    ML_mod(1,:)=0; ML_mod(1,1)=1;
end
uFCT = ML_mod \ rhs;

end