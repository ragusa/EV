function D = compute_high_order_diffusion_matrix(u_older,u_old,dt_old,...
    viscL,mesh,phys,quadrature,ev,dof_handler,high_order_scheme)

% compute entropy viscosity
if (high_order_scheme == 2) % Entropy viscosity
    viscE = compute_entropy_viscosity(...
        u_older,u_old,dt_old,mesh,phys,quadrature,ev,dof_handler);
else % No viscosity
    viscE = zeros(mesh.n_cell,1);
end

% compute high-order viscosity
viscH = min(viscE,viscL);

% compute diffusion matrix
D = compute_diffusion_matrix(viscH,mesh.dx,dof_handler);

end
