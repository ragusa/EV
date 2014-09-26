function [MC,ML,K,D,b]=build_matrices(len,nel,omega,sigma,src,speed,inc,impose_strongly)
% builds the mass matrices, steady state system matrix, artificial
% dissipation matrix, and steady state rhs
%
% MC = consistent mass matrix
% ML = lumped mass matrix
% K  = Kuzmin's K matrix
% D  = Kuzmin's D matrix
% b  = steady-state rhs

% integral bi bj
M_cell=[2 1;1 2]/3/speed;
% integral omega * bi * dbj/dx
if (impose_strongly)
    K_cell= omega*[-1 1;-1 1]/2;
else
    K_cell=-omega*[-1 1;-1 1]/2; K_cell = K_cell';
end
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
AH = spalloc(n,n,nnz);
D  = spalloc(n,n,nnz);
b  = zeros(n,1);

% compute sigma for each cell
sigma_cell = zeros(nel,1);
% center x-value for each cell
x_center = linspace(0.5*h,len-0.5*h,nel);
for iel = 1:nel
    sigma_cell(iel) = sigma(x_center(iel));
end

% assemble consistent mass matrix and steady-state rhs
for iel=1:nel
    MC(g(iel,:),g(iel,:)) = MC(g(iel,:),g(iel,:)) + M_cell*jac;
    AH(g(iel,:),g(iel,:)) = AH(g(iel,:),g(iel,:)) + K_cell + sigma_cell(iel)*M_cell*jac;
    b(g(iel,:))           = b(g(iel,:))           + b_cell*src*jac;
end
if (~impose_strongly)
    AH(end,end) = AH(end,end) + omega;
end

% kuzmin writes K on the rhs
K = -AH;

% lump mass
ML=(sum(MC));

% compute D
for i=1:n
    for j=1:n
        D(i,j)=max([0,AH(i,j),AH(j,i)]);
    end
    D(i,i)=-(sum(D(i,1:i-1))+sum(D(i,i+1:n)));
end

if (~impose_strongly)
    b(1)=b(1)+omega*inc;
end

return
end
