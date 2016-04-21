function [f_min,f_max] = compute_min_max_per_dof(f,t,n_dof,mesh,zq,s)

% extract mesh
x      = mesh.x;
n_cell = mesh.n_cell;

% compute range of min/max operations
range = ceil(s/mesh.dx(1));

% compute min and max of function in the support of i
f_min  = 1e15*ones(n_dof,1);
f_max  = zeros(n_dof,1);
for iel = 1:n_cell
    % compute quadrature point positions
    xq = get_quadrature_point_positions(x,iel,zq);
    
    % compute function at each quadrature point
    f_cell = f(xq,t);
    
    % compute max and min on cell
    f_cell_min = min(f_cell);
    f_cell_max = max(f_cell);
    
    % update max/min for each dof on cell
    iL = max(iel-range+1,1);
    iR = min(iel+range,n_dof);
    f_min(iL:iR) = min(f_min(iL:iR),f_cell_min);
    f_max(iL:iR) = max(f_max(iL:iR),f_cell_max);
end

% for ranges that extend beyond first node, set the minimum to zero
% because both sigma and source should be considered zero when imposing
% Dirichlet BC
for i = 2:range
    f_min(i) = 0.0;
end

end