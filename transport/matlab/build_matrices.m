%%
% Description:
%   Builds the mass matrices, inviscid and low-order steady-state matrices,
%   low-order viscosity and corresponding artificial diffusion matrix.
%
% MC = consistent mass matrix
% ML = lumped mass matrix
% A  = steady-state matrix
% D  = artificial diffusion matrix

%%
function [MC,ML,A,AL,D,viscL] = build_matrices(...
    nq,zq,wq,v,dvdz,Jac,x,dx,nel,n,connectivity,...
    mu,speed,sigma,low_order_scheme,modify_for_weak_DirichletBC)

%--------------------------------------------------------------------------
% Compute mass matrices and inviscid steady-state matrix
%--------------------------------------------------------------------------

% degrees of freedom per cell
dofs_per_cell = 2;

% intialize matrices
nnz = 3*n; % max number of nonzero entries
MC = spalloc(n,n,nnz);
A  = spalloc(n,n,nnz);
D  = spalloc(n,n,nnz);

% assemble consistent mass matrix and inviscid steady-state system matrix
for iel = 1:nel
    % get local dof indices
    ii = connectivity(iel,:);
    
    % compute quadrature point positions
    xq = get_quadrature_point_positions(x,iel,zq);
    
    % assemble local matrices
    M_cell = zeros(dofs_per_cell,dofs_per_cell);
    A_cell = zeros(dofs_per_cell,dofs_per_cell);
    for q = 1:nq
        M_cell = M_cell + v(:,q) * v(:,q)' * wq(q) * Jac(iel);
        A_cell = A_cell + mu * v(:,q) * dvdz(:,q)' * wq(q) +...
            sigma(xq(q)) * v(:,q) * v(:,q)' * wq(q) * Jac(iel);
    end
    
    % add local contributions to global system
    MC(ii,ii) = MC(ii,ii) + M_cell;
    A(ii,ii)  = A(ii,ii)  + A_cell;
end

% subtract inflow boundary face term if weakly imposing Dirichlet BC
if modify_for_weak_DirichletBC
    A(1,1) = A(1,1) + 0.5 * mu;
end

% multiply steady-state matrix by speed
A = A*speed;

% computed lumped mass matrix
ML = diag(sum(MC));

%--------------------------------------------------------------------------
% Compute low-order viscosity
%--------------------------------------------------------------------------

% compute sum of local viscous bilinear forms
visc_bilin_cell = [1 -1; -1 1];
visc_bilin = zeros(n,n);
for iel = 1:nel
    % get local dof indices
    ii = connectivity(iel,:);
    
    % add local contributions to global system
    visc_bilin(ii,ii) = visc_bilin(ii,ii) + dx(iel)*visc_bilin_cell;
end

% compute low-order viscosity
viscL = zeros(nel,1);
for iel = 1:nel
    viscL(iel) = 0; % initialize to zero for max()
    for i = 1:dofs_per_cell
        for j = 1:dofs_per_cell
            if (i ~= j)
                ii = connectivity(iel,i);
                jj = connectivity(iel,j);
                viscL(iel) = max(viscL(iel),-max(0,A(ii,jj))/visc_bilin(ii,jj));
            end
        end
    end
end

%--------------------------------------------------------------------------
% Compute low-order artificial diffusion matrix and low-order steady-state
% matrix
%--------------------------------------------------------------------------

% compute low-order artificial diffusion matrix
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
            % get local dof indices
            ii = connectivity(iel,:);
    
            % add local contributions to global system
            D(ii,ii) = D(ii,ii) + viscL(iel) * dx(iel)*visc_bilin_cell;
        end
    otherwise
        error('Invalid low order scheme chosen');
end

% compute low-order steady-state matrix
AL = A + D;

return
end