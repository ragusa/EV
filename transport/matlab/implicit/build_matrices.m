function [MC,ML,K,D,b]=build_matrices(len,nel,omega,sigma,src)
% builds the mass matrices, steady state system matrix, artificial
% dissipation matrix, and steady state rhs
%
% MC = consistent mass matrix
% ML = lumped mass matrix
% K  = Kuzmin's K matrix
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
    K(g(iel,:),g(iel,:))  = K(g(iel,:),g(iel,:))  + K_cell;
    b(g(iel,:))           = b(g(iel,:))           + b_cell*src*jac;
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
