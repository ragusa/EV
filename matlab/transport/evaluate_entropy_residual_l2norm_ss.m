function l2norm = evaluate_entropy_residual_l2norm_ss(...
    u,mesh,phys,quadrature,ev,dof_handler)

% unpack quadrature
zq   = quadrature.zq;
wq   = quadrature.wq;
v    = quadrature.v;
dvdz = quadrature.dvdz;
Jac  = quadrature.Jac;

% unpack physics
mu     = phys.mu;
sigma  = phys.sigma;
source = phys.source;

% unpack mesh
x   = mesh.x;
nel = mesh.n_cell;

% unpack dof handler
n_dof = dof_handler.n_dof;
g     = dof_handler.connectivity;

% unpack entropy viscosity
entropy       = ev.entropy;
entropy_deriv = ev.entropy_deriv;

% compute integral
integral = 0.0;
for iel = 1:nel
    % evaluate new and old solution at quadrature points
    u_local_nodes = u(g(iel,:));
    u_local = v' * u_local_nodes;
    dudx_local = dvdz' * u(g(iel,:)) / Jac(iel);
    % evaluate entropy and derivative
    E    = entropy(u_local);
    dEdu = entropy_deriv(u_local);
    dEdx = dEdu .* dudx_local;
    % compute quadrature point positions
    xq = get_quadrature_point_positions(x,iel,zq);
    % evaluate cross section and source at quadrature points
    sigma_q  = sigma(xq);
    source_q = source(xq);
    
    % compute maximum entropy residual at quadrature points
    R_q = mu*dEdx + dEdu.*(sigma_q.*u_local - source_q);

    % add to integral
    integral = integral + wq' * R_q.^2 * Jac(iel);
end

% compute square root to finish l2 norm
l2norm = sqrt(integral);

