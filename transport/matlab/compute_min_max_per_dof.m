function [f_min,f_max] = compute_min_max_per_dof(f,t,n_dof,mesh,zq)

% extract mesh
x      = mesh.x;
n_cell = mesh.n_cell;

% compute min and max sigma and source in the support of i
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
    f_min(iel:iel+1) = min(f_min(iel:iel+1),f_cell_min);
    f_max(iel:iel+1) = min(f_max(iel:iel+1),f_cell_max);
end

end