function [D,viscE] = compute_high_order_diffusion_matrix(u_older,u_old,dt_old,...
    viscL,mesh,phys,quadrature,ev,dof_handler,high_order_scheme)

% compute entropy viscosity
if (high_order_scheme == 0)
    % no viscosity
    viscE = zeros(mesh.n_cell,1);
else
    % entropy viscosity
    viscE = compute_entropy_viscosity_alternate(...
        u_older,u_old,dt_old,mesh,phys,quadrature,ev,dof_handler,...
        high_order_scheme);
end

% compute high-order viscosity
viscH = min(viscE,viscL);

% compute diffusion matrix
D = compute_diffusion_matrix(viscH,mesh.dx,dof_handler);

end
