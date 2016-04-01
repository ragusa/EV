function [uH,AH,DH] = compute_high_order_solution_theta(u_old,dt,MC,...
    A,AH,b_old,b_new,viscL,quadrature,mesh,dof_handler,phys,ev,numerics_opts,...
    nonlin_opts)

% extract options
modify_for_strong_DirichletBC = numerics_opts.modify_for_strong_DirichletBC;
high_order_scheme             = numerics_opts.high_order_scheme;
theta                         = numerics_opts.theta;
tol                  = nonlin_opts.nonlin_tol;
max_iter             = nonlin_opts.max_iter;
relaxation_parameter = nonlin_opts.relax;

% initialize solution iterate
uH = u_old;

% compute old ss residual
ss_res = b_old - AH*u_old;

% compute system rhs
system_rhs = MC*u_old + (1-theta)*dt*ss_res ...
    + theta*dt*b_new;
if (modify_for_strong_DirichletBC)
    system_rhs(1) = phys.inc;
end

% begin nonlinear iteration
for iter = 1:max_iter
    % compute high-order diffusion matrix
    DH = compute_high_order_diffusion_matrix(u_old,...
        uH,dt,viscL,mesh,phys,quadrature,ev,...
        dof_handler,high_order_scheme);
    
    % compute high-order ss matrix
    AH = A + DH;
    
    % compute system matrix
    system_matrix = MC + theta*dt*AH;
    if (modify_for_strong_DirichletBC)
        system_matrix(1,:)=0; system_matrix(1,1)=1;
    end
    
    % compute residual
    res = system_rhs - system_matrix*uH;

    % test convergence of previous iteration
    converged = test_convergence(res,tol,iter);
    if (converged)
        fprintf('\tConverged at iteration %i\n',iter);
        break;
    end
    
    % compute change in solution iterate
    du = system_matrix \ system_rhs - uH;
    
    % update solution
    uH = uH + relaxation_parameter * du;
end

% report if the solution did not converge
if (~converged)
    error('High-order Solution did not converge');
end

end