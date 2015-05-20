function converged = test_convergence(res,tol,iter)

resnorm = norm(res,2);
fprintf('\tIteration %i: residual error = %e\n',iter,resnorm);
converged = resnorm < tol;


end