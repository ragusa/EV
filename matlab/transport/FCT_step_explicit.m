function [uFCT,Wminus,Wplus] = FCT_step_explicit(u_old,uH,dt,ML,MC,AL,DH,DL,b,...
    sigma_min,sigma_max,source_min,source_max,mesh,phys,fct_opts,opts)

% unpack options
prelimit = fct_opts.prelimit;
modify_for_strong_DirichletBC = opts.modify_for_strong_DirichletBC;
inc = phys.inc;

% size of system
n_dof = length(u_old);

% compute flux corrections
F = flux_correction_matrix(u_old,uH,dt,DH,DL,MC,0.0);

% compute low-order solution
uL = compute_low_order_solution_theta(u_old,AL,...
  ML,b,dt,0,inc,modify_for_strong_DirichletBC);

% prelimit flux corrections if user specified
if (prelimit)
    % prelimit fluxes
    F = prelimit_fluxes(F,uL);
end
        
% compute limited antidiffusion source
[flim,Wminus,Wplus] = compute_limited_flux_sums(...
    u_old,u_old,uL,dt,ML,AL,b,F,...
    sigma_min,sigma_max,source_min,source_max,mesh,0.0,n_dof,...
    phys,fct_opts,opts);

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
