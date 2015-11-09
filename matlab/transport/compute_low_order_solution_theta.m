function uL = compute_low_order_solution_theta(u_old,AL,ML,...
    b,dt,theta,inc,modify_for_strong_DirichletBC)

% compute system matrix and rhs
system_matrix = ML + dt*theta*AL;
system_rhs = (ML - (1-theta)*dt*AL)*u_old + dt*b;
if modify_for_strong_DirichletBC
    system_matrix(1,:)=0; system_matrix(1,1)=1;
    system_rhs(1) = inc;
end

% compute solution
uL = system_matrix \ system_rhs;

end