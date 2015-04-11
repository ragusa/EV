%%
% Description:
%   Assembles the steady-state right hand side vector, which may be a
%   function of time.
%
% Input:
%   [name]
%
% Ouptut:
%   [name]
%

%%
function b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n,connectivity,...
    source,mu,speed,t,inc,modify_for_weak_DirichletBC)

% degrees of freedom per cell
dofs_per_cell = 2;

% intialize
b  = zeros(n,1);

% assemble consistent mass matrix and steady-state rhs
for iel = 1:nel
    % get local dof indices
    ii = connectivity(iel,:);
    
    % compute quadrature point positions
    xq = get_quadrature_point_positions(x,iel,zq);
    
    % assemble local rhs
    b_cell = zeros(dofs_per_cell,1);
    for q = 1:nq
        b_cell = b_cell + source(xq(q),t) * v(:,q) * wq(q) * Jac(iel);
    end
    
    % add local contributions to global system
    b(ii) = b(ii) + b_cell;
end

% if weakly imposing Dirichlet BC, add inflow boundary term
if modify_for_weak_DirichletBC
    b(1) = b(1) + 0.5 * mu * inc;
end

% multiply by speed
b = b*speed;

end