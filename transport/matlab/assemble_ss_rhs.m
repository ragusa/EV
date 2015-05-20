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
function b = assemble_ss_rhs(t,quadrature,mesh,dof_handler,phys,...
    modify_for_weak_DirichletBC)

% unpack quadrature
nq   = quadrature.nq;
zq   = quadrature.zq;
wq   = quadrature.wq;
v    = quadrature.v;
Jac  = quadrature.Jac;

% unpack physics
mu     = phys.mu;
source = phys.source;
speed  = phys.speed;
inc    = phys.inc;

% unpack mesh
x   = mesh.x;
nel = mesh.n_cell;

% unpack dof handler
n = dof_handler.n_dof;
g = dof_handler.connectivity;
dofs_per_cell = dof_handler.dofs_per_cell;



% intialize
b  = zeros(n,1);

% assemble consistent mass matrix and steady-state rhs
for iel = 1:nel
    % get local dof indices
    ii = g(iel,:);
    
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