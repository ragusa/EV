function F = compute_flux_correction_matrix_implicit(u_old,uH,dt,D,MC)
% computes the flux correction matrix
%
% F     = flux correction matrix
%
% u_old = old solution
% uH    = high-order solution
% dt    = time step size
% D     = Kuzmin's artificial dissipation matrix

% size of system
n = length(u_old);

% compute flux correction matrix
F = zeros(n,n);
for i = 1:n
    for j = 1:n
        F(i,j) = -MC(i,j)*(uH(j)-u_old(j) - (uH(i)-u_old(i))) ...
            -dt*D(i,j)*(uH(j)-uH(i));
    end
end

return
end