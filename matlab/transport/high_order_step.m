function [uH,DH] = high_order_step(u_older,u_old,dt_old,dt,A,b,MC,theta,...
    viscL,mesh,phys,quadrature,dof_handler,ev,...
    high_order_scheme,modify_for_strong_DirichletBC)

% compute high-order diffusion matrix
DH = compute_high_order_diffusion_matrix(u_older,...
    u_old,dt_old,viscL,mesh,phys,quadrature,ev,...
    dof_handler,high_order_scheme);

% compute high-order matrices
AH = A + DH;
system_matrix = MC + dt*theta*AH;
if modify_for_strong_DirichletBC
    system_matrix(1,:)=0; system_matrix(1,1)=1;
end

% compute high-order solution
system_rhs = (MC - (1-theta)*dt*AH)*u_old + dt*b;
if modify_for_strong_DirichletBC
    system_rhs(1) = phys.inc;
end
uH = system_matrix \ system_rhs;

end