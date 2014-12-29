function [MC,ML,A,AL,D,b,viscL] = ...
    build_matrices(len,nel,omega,speed,sigma,source,...
    low_order_scheme,high_order_scheme,periodic_BC)
% builds the mass matrices, steady state system matrix, artificial
% dissipation matrix, and steady state rhs
%
% MC = consistent mass matrix
% ML = lumped mass matrix
% K  = Kuzmin's K matrix
% D  = Kuzmin's D matrix
% b  = steady-state rhs

% local bilinear forms (cells are uniform)
M_cell = [2 1;1 2]/3;               % v(i)*v(j)
R_cell = speed*M_cell;              % c*v(i)*v(j)
K_cell = speed*omega*[-1 1;-1 1]/2; % c*omega*v(i)*dv(j)/dx
b_cell = speed*[1;1];               % c*v(i)

h = len/nel; % cell size
jac = h/2;   % Jacobian of transformation to reference cell
g = [linspace(1,nel,nel)' linspace(2,nel+1,nel)']; % connectivity array
if periodic_BC
    g(end,:) = [g(end,1) 1];
end
dofs_per_cell = 2; % degrees of freedom per cell

% intialize matrices
n = nel+1; % system size
if periodic_BC, n=nel; end
nnz = 3*n; % max number of nonzero entries
MC = spalloc(n,n,nnz);
A  = spalloc(n,n,nnz);
D  = spalloc(n,n,nnz);
b  = zeros(n,1);

% assemble consistent mass matrix and steady-state rhs
for iel = 1:nel
    MC(g(iel,:),g(iel,:)) = MC(g(iel,:),g(iel,:)) + M_cell*jac;
    A(g(iel,:),g(iel,:))  = A(g(iel,:),g(iel,:))  + K_cell + sigma(iel)*R_cell*jac;
    b(g(iel,:))           = b(g(iel,:))           + b_cell*source(iel)*jac;
end

% lump mass matrix and reaction matrix
ML = diag(sum(MC));

% compute low-order viscosity if needed
if low_order_scheme == 2 || high_order_scheme == 2
    % compute sum of local viscous bilinear forms
    visc_bilin_cell = h*[1 -1; -1 1];
    visc_bilin = zeros(n,n);
    for iel = 1:nel
        visc_bilin(g(iel,:),g(iel,:)) = visc_bilin(g(iel,:),g(iel,:))...
            + visc_bilin_cell;
    end
    
    % compute low-order viscosity
    viscL = zeros(nel,1);
    for iel = 1:nel
        viscL(iel) = 0; % initialize to zero for max()
        for i = 1:dofs_per_cell
            for j = 1:dofs_per_cell
                if (i ~= j)
                    viscL(iel) = max(viscL(iel),-abs(A(i,j))/visc_bilin(i,j));
                end
            end
        end
    end
else
    % pass an empty container if low-order viscosity isn't needed
    viscL = [];
end

% compute discrete upwinding operator D for low-order advection matrix
switch low_order_scheme
    case 1 % algebraic
        for i = 1:n
            for j = 1:n
                D(i,j) = -max([0,A(i,j),A(j,i)]);
            end
            D(i,i) = -(sum(D(i,1:i-1))+sum(D(i,i+1:n)));
        end
    case 2 % graph-theoretic
        for iel = 1:nel
            D(g(iel,:),g(iel,:)) = D(g(iel,:),g(iel,:)) + viscL(iel)*visc_bilin_cell;
        end
    otherwise
        error('Invalid low order scheme chosen');
end

% compute low-order and high-order system matrices
AL = A + D;

return
end
