function value = myfunc(L)

% number of cells, dofs, edges
n_cells = 32;
n_dofs = n_cells + 1;
n_edges = n_dofs - 1;

% mesh
x = linspace(0.0,1.0,n_cells+1)'; % positions of dofs/cell vertices
dx = diff(x);

% quadrature points and test functions
nq = 4;
[zq,wq] = get_GL_quadrature(nq);
[v,dvdz] = get_lagrange_basis(zq,2); 
Jac = 0.5*dx;

% get low-order solution
x_uL = importdata('uL_ss.txt');
uL = x_uL(:,2);

% get fluxes
flux_matrix = importdata('flux.txt');
flux = zeros(n_edges,1);
for iedge = 1:n_edges
    flux(iedge) = flux_matrix(iedge,iedge+1);
end

% compute u^{n+1}
u = compute_uFCT(uL,L,flux);

value = 0.0;
for iel = 1:n_cells
    % compute entropy at quadrature points
    u_left = u(iel);
    u_right = u(iel+1);
    uq = v(1,:)*u_left + v(2,:)*u_right;
    entropy = 0.5*uq.*uq;
    
    % compute integrals
    for i = 1:2
        value = value + (v(i,:).*entropy)*(Jac(iel)*wq);
    end
end

return;
end