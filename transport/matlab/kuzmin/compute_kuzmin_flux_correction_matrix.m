function F = compute_kuzmin_flux_correction_matrix(u_old,uH,dt,D,MC,use_explicit)
% computes Kuzmin's flux correction matrix for either explicit or implicit
% Euler
%
% F     = flux correction matrix
%
% u_old = old solution
% uH    = high-order solution
% dt    = time step size
% D     = Kuzmin's artificial dissipation matrix
% use_explicit = boolean for whether explicit Euler is used (as opposed to
%                implicit Euler)

% size of system
n = length(u_old);

% compute flux correction matrix
F = zeros(n,n);
if (use_explicit) % explicit Euler
    for i = 1:n
        for j = 1:n
            F(i,j) = -MC(i,j)*(uH(j)-u_old(j) - (uH(i)-u_old(i))) ...
                -dt*D(i,j)*(u_old(j)-u_old(i));
        end
    end
else              % implicit Euler
    for i = 1:n
        for j = 1:n
            F(i,j) = -MC(i,j)*(uH(j)-u_old(j) - (uH(i)-u_old(i))) ...
                -dt*D(i,j)*(uH(j)-uH(i));
        end
    end
end

return
end