function D = compute_diffusion_matrix(visc,dx,dof_handler)

% unpack
n   = dof_handler.n_dof;
nel = dof_handler.n_cell;
g   = dof_handler.connectivity;

% initialize
D  = spalloc(n,n,3*n);

% cell viscous bilinear form
visc_bilin_cell = [1 -1; -1 1];

% loop over cells
for iel = 1:nel
    % get local dof indices
    ii = g(iel,:);
    
    % add local contributions to global system
    D(ii,ii) = D(ii,ii) + visc(iel) * dx(iel)*visc_bilin_cell;
end

end