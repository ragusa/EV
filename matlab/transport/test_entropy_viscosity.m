% input parameters
nel = 32;
x = linspace(0,10,nel+1)';
u_new = sin(2*x);
u_old = sin(x);
mu = 0.3;
sigma = 2*ones(nel,1);
source = 3*ones(nel,1);
dt = 0.01;
nq = 3;
cE = 1.0;
cJ = 0.0;
entropy = @(u) 0.5*u.*u;
entropy_deriv = @(u) u;
periodic_BC = false;

[zq,wq]  = get_GL_quadrature(nq);
[v,dvdz] = get_lagrange_basis(zq,2);
dx_cell = diff(x);
Jac = 0.5*dx_cell;

% compute entropy viscosity
viscE = compute_entropy_viscosity(...
    u_old,u_new,x,mu,sigma,source,dt,...
    v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,periodic_BC);

% put entropy viscosity data in deal.ii output format
x_sav  = zeros(nel*3,1);
ev_sav = zeros(nel*3,1);
for i = 1:nel
    xL = x(i);
    xR = x(i+1);
    xC = 0.5*(xL+xR);
    x_sav((i-1)*3+1) = xL;
    x_sav((i-1)*3+2) = xC;
    x_sav((i-1)*3+3) = xR;
    ev_sav((i-1)*3+1) = viscE(i);
    ev_sav((i-1)*3+2) = viscE(i);
    ev_sav((i-1)*3+3) = viscE(i);
end

% plot entropy viscosity
plot(x_sav,ev_sav,'-x');

% save entropy viscosity
dlmwrite('output/entropy_viscosity_MATLAB.txt',[x_sav,ev_sav],' ');

% save solution
xq = zeros(nq*nel,1);
u_new_q = zeros(nq*nel,1);
u_old_q = zeros(nq*nel,1);
for iel=1:nel
    for q = 1:nq
        xq(nq*(iel-1)+q) = x(iel) + (zq(q)+1)/2*dx_cell(iel);
        u_new_q(nq*(iel-1)+q) = v(1,q)*u_new(iel) + v(2,q)*u_new(iel+1);
        u_old_q(nq*(iel-1)+q) = v(1,q)*u_old(iel) + v(2,q)*u_old(iel+1);
    end
end
dlmwrite('output/solutions_MATLAB.txt',[xq,u_new_q,u_old_q],' ');