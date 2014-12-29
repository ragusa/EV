function D = compute_high_order_diffusion_matrix(viscE,viscL,h,periodic_BC)

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
D_cell = h*[1 -1; -1 1]; % local bilinear form
D  = spalloc(n,n,nnz);
for iel = 1:nel
    viscH = min(viscL(iel),viscE(iel));
    D(g(iel,:),g(iel,:)) = D(g(iel,:),g(iel,:)) + viscH*D_cell;
end

end
