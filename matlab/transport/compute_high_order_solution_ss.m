function [uH,DH,viscE] = compute_high_order_solution_ss(A,b,viscL,mesh,...
    phys,quadrature,ev,dof_handler,opts,nonlin_opts,out_opts)

% unpack options
high_order_scheme = opts.high_order_scheme;
modify_for_strong_DirichletBC = opts.modify_for_strong_DirichletBC;
max_iter = nonlin_opts.max_iter;
nonlin_tol = nonlin_opts.nonlin_tol;
relaxation_parameter = nonlin_opts.relax;
plot_iterations = out_opts.plot_EV_iteration;

% initialize solution iterate
uH = zeros(dof_handler.n_dof,1);

% begin nonlinear iteration
converged = false;
for iter = 1:max_iter
    % compute high-order diffusion matrix
    [DH,viscE] = compute_high_order_diffusion_matrix(uH,...
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

    % enforce Dirichlet BC; this is not necessary - it should already be done
    if (modify_for_strong_DirichletBC)
      uH(1) = phys.inc;
    end

    % plot iteration if requested
    if (plot_iterations)
        % new figure
        figure(1);

        % plot solution
        subplot(2,1,1);
        plot(mesh.x, uH);
        leg_string = sprintf('High-order, iteration %i',iter);
        legend(leg_string);

        % plot viscosity
        subplot(2,1,2);
        semilogy(mesh.x_center,viscE);
        leg_string = sprintf('Entropy viscosity, iteration %i',iter);
        legend(leg_string,'Location','NorthWest');

        % pause
        if (strcmp(out_opts.pause_type,'wait'))
            waitforbuttonpress;
        else
            pause(out_opts.pausetime);
        end
        
    end
end

% report if the solution did not converge
if (~converged)
    error('High-order Solution did not converge');
end

end
