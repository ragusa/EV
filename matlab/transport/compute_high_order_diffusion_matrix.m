function [D,viscE] = compute_high_order_diffusion_matrix(u_older,u_old,dt_old,...
    viscL,mesh,phys,quadrature,ev,dof_handler,high_order_scheme)

% compute entropy viscosity
if (high_order_scheme == 2) % original entropy viscosity
    viscE = compute_entropy_viscosity(...
        u_older,u_old,dt_old,mesh,phys,quadrature,ev,dof_handler);
elseif (high_order_scheme == 3) % alternate entropy viscosity
    viscE = compute_entropy_viscosity_alternate(...
        u_older,u_old,dt_old,mesh,phys,quadrature,ev,dof_handler);
elseif (high_order_scheme == 4) % alternate entropy viscosity 2
    viscE = compute_entropy_viscosity_alternate2(...
        u_older,u_old,dt_old,mesh,phys,quadrature,ev,dof_handler);
else % No viscosity
    viscE = zeros(mesh.n_cell,1);
end

% compute high-order viscosity
viscH = min(viscE,viscL);

% compute diffusion matrix
D = compute_diffusion_matrix(viscH,mesh.dx,dof_handler);

end
