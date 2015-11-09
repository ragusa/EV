function F = flux_correction_matrix(u_old,uH,dt,DH,DL,MC,theta)
% computes Kuzmin's flux correction matrix for either explicit or implicit
% Euler
%
% F     = flux correction matrix
%
% u_old = old solution
% uH    = high-order solution
% dt    = time step size
% D     = Kuzmin's artificial dissipation matrix
% MC    = consistent mass matrix
% theta = time-discretization parameter for theta-scheme

% size of system
n = length(u_old);

% change in solution
du = uH - u_old;

% compute flux correction matrix
F = zeros(n,n);
for i = 1:n
    for j = 1:n
        F(i,j) = -MC(i,j)*(du(j) - du(i))/dt ...
            + (DL(i,j)-DH(i,j))*(theta*(uH(j)-uH(i)) + (1-theta)*(u_old(j)-u_old(i)));
    end
end

return
end