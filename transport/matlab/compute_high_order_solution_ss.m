function [uH,DH] = compute_high_order_solution_ss(A,b,viscL,mesh,...
    phys,quadrature,ev,dof_handler,high_order_scheme,max_iter,...
    nonlin_tol,relaxation_parameter,modify_for_strong_DirichletBC)

% initialize solution iterate
uH = zeros(dof_handler.n_dof,1);

% begin nonlinear iteration
for iter = 1:max_iter
    % compute high-order diffusion matrix
    DH = compute_high_order_diffusion_matrix(uH,...
        uH,1.0,viscL,mesh,phys,quadrature,ev,...
        dof_handler,high_order_scheme);
    
    % compute system matrix and rhs
    system_matrix = A + DH;
    system_rhs = b;
    if (modify_for_strong_DirichletBC)
        system_matrix(1,:)=0; system_matrix(1,1)=1;
        system_rhs(1) = phys.inc;
    end
    
    % compute residual
    res = system_rhs - system_matrix*uH;

    % test convergence of previous iteration
    converged = test_convergence(res,nonlin_tol,iter);
    if (converged)
        fprintf('\tConverged at iteration %i\n',iter);
        break;
    end
    
    % update solution
    uH = uH + relaxation_parameter * (system_matrix \ res);
end

% report if the solution did not converge
if (~converged)
    error('\t\tHigh-order Solution did not converge\n');
end

end
