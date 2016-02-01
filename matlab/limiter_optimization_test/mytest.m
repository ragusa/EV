clear; clc;

% number of cells, dofs, edges
n_cells = 32;
n_dofs = n_cells + 1;
n_edges = n_dofs - 1;

% mesh
x = linspace(0.0,1.0,n_cells+1)'; % positions of dofs/cell vertices

% get low-order steady-state matrix
AL = importdata('AL.txt');
AL(1,:) = 0; AL(:,1) = 0; AL(1,1) = 1.0;

% get fluxes
flux_matrix = importdata('flux.txt');
flux = zeros(n_edges,1);
for iedge = 1:n_edges
    flux(iedge) = flux_matrix(iedge,iedge+1);
end

% compute optimal limiting coefficients
A = zeros(n_edges*2,n_edges);
b = zeros(n_edges*2,1);
for i = 1:n_edges
    A(i*2-1,i) = 1.0;
    A(i*2,i) = -1.0;
    b(i*2-1) = 1.0;
    b(i*2) = 0.0;
end
guess = zeros(n_edges,1);
L = fmincon(@myfunc,guess,A,b);

%% Plot
figure; clf; hold on;

% plot low-order solution
x_uL = importdata('uL_ss.txt');
uL = x_uL(:,2);
plot(x,uL,'r-x');

% plot high-order solution
L_matrix = ones(n_dofs,n_dofs);
limited_fluxes = sum(flux_matrix.*L_matrix,2); limited_fluxes(1) = 0.0;
du = AL\limited_fluxes;
uH = uL + du;
plot(x,uH,'b-+');

% plot FCT solution
uFCT = compute_uFCT(uL,L,flux);
plot(x,uFCT,'g-o');

% plot exact solution
x_u = importdata('uexact.txt');
plot(x_u(:,1),x_u(:,2),'k-');
legend('Low','High','FCT','Exact');
