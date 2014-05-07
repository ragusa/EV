function [MC,K,ML,D,b]=build_matrices(len,nel,omega,sigma,src,inc,impose_weakly)
% builds the mass matrices, steady state system matrix, artificial
% dissipation matrix, and steady state rhs
%
% MC = consistent mass matrix
% K  = Kuzmin's K matrix
% ML = lumped mass matrix
% D  = Kuzmin's D matrix
% b  = steady-state rhs

% integral bi bj
M_cell=[2 1;1 2]/3;
% integral omega * bi * dbj/dx
K_cell=omega*[-1 1;-1 1]/2;
% integral -omega * dbi/dx bj
K_cell_weak = -K_cell';
% integral bi
b_cell=[1;1];

h=len/nel;
jac=h/2;
% connectivity array CFEM
g=[linspace(1,nel,nel)' linspace(2,nel+1,nel)'];

% intialize matrices
n = nel+1; % system size
nnz = 3*n; % max number of nonzero entries
MC = spalloc(n,n,nnz);
K  = spalloc(n,n,nnz);
D  = spalloc(n,n,nnz);
b  = zeros(n,1);

% assemble consistent mass matrix and steady-state rhs
for iel=1:nel
    MC(g(iel,:),g(iel,:)) = MC(g(iel,:),g(iel,:)) + M_cell*jac;
    b(g(iel,:))           =  b(g(iel,:))          + b_cell*src*jac;
end

% assemble K matrix and modify rhs if necessary for Dirichlet BC
if (impose_weakly) % build matrix with weakly imposed Dirichlet BC
    % add interior contributions
    for iel = 1:nel
        K(g(iel,:),g(iel,:)) = K(g(iel,:),g(iel,:)) + K_cell_weak;
    end
    % add outflow contribution
    K(n,n) = K(n,n) + omega;
    % add inflow contribution to rhs
    b(1) = b(1) + omega*inc;
else % build matrix without doing any integration by parts
    for iel = 1:nel
        K(g(iel,:),g(iel,:)) = K(g(iel,:),g(iel,:)) + K_cell;
    end
end

% kuzmin writes K on the rhs
K=-K-sigma*MC;

% lump mass
ML=(sum(MC));

% compute D
for i=1:n
    for j=1:n
        D(i,j)=max([0,-K(i,j),-K(j,i)]);
    end
    D(i,i)=-(sum(D(i,1:i-1))+sum(D(i,i+1:n)));
end

return
end
