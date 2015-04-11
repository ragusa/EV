function D = compute_high_order_diffusion_matrix(viscE,viscL,dx,periodic_BC)

nel = length(viscE); % number of elements
n = nel+1;           % number of dofs
if periodic_BC
    n = nel;
end
nnz = 3*n;           % max number of nonzero entries
g = [linspace(1,nel,nel)' linspace(2,nel+1,nel)']; % connectivity array
if periodic_BC
    g(end,:) = [ g(end,1) 1];
end

% build matrix
D  = spalloc(n,n,nnz);
for iel = 1:nel
    % get local dof indices
    ii = g(iel,:);
    
    % compute cell matrix
    viscH = min(viscL(iel),viscE(iel));
    D_cell = dx(iel)*viscH*[1 -1; -1 1];    
    
    % add local contribution to global system
    D(ii,ii) = D(ii,ii) + D_cell;
end

end
