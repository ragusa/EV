function [MC,K,ML,D,b]=build_matrices(len,nel,omega,sigma,src)

% integral bi bj
m=[2 1;1 2]/3;
% integral bi dbj/dx
gr=[-1 1;-1 1]/2;
% integral bi
f=[1;1];

h=len/nel;
jac=h/2;
% connectivity array CFEM
g=[linspace(1,nel,nel)' linspace(2,nel+1,nel)'];

% intialize matrices
n=nel+1;
nnz_=3*n;
MC=spalloc(n,n,nnz_);
K=MC;
D=MC;
b=zeros(n,1);

for iel=1:nel
    MC(g(iel,:),g(iel,:)) = MC(g(iel,:),g(iel,:)) + jac*m;
    K(g(iel,:),g(iel,:))  = K(g(iel,:),g(iel,:))  + gr;
    b(g(iel,:)) = b(g(iel,:)) +jac*f*src;
end

% kuzmin write K on the rhs
K=-omega*K-sigma*MC;

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
