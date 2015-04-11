%%
% Description:
%   Computes the quadrature point positions in a cell.
%
% Input:
%   [name]  [size]   [description]
%   x       n        cell vertices for all n cells
%   i       1        cell index (assumed to start from 1)
%   zq      nq       quadrature point positions in reference cell (-1,1)
%
% Output:
%   [name]  [size]  [description]
%   xq      nq      quadrature point positions in cell i

%%
function xq = get_quadrature_point_positions(x,i,zq)

% compute Jacobian of cell
Jac = 0.5*(x(i+1)-x(i));

% compute quadrature point positions
xq = x(i) + Jac*(1 + zq);

end