function F = compute_flux_correction_matrix(u_old,dt,D,B,K,b)
% computes the flux correction matrix:
%    sum(F,2) = -dt*D*u_old + dt*B*(K*u_old + b), where sum(F,2) is row-sum
% This is decomposed as
%    F_ij = dt*D_ij*(u_old,j - u_old,i) + dt*(B_ji*G_i - B_ij*G_j)
%
% F     = flux correction matrix
% u_old = old solution
% dt    = time step size
% D     = Kuzmin's artificial dissipation matrix
% B     = Guermond's B matrix: B = (ML-MC)*ML^-1
% K     = Kuzmin's K matrix; negative of high-order steady-state matrix
% b     = steady-state rhs vector

% size of system
n = length(u_old);

% compute the G vector in Guermond's paper:
% MC*du/dt = -G          (Guermond)
% MC*du/dt = K*u_old + b (Kuzmin)
%    => G = -K*u_old - b
G = -K*u_old - b;

% compute flux correction matrix
F = zeros(n,n);
for i = 1:n
    for j = 1:n
        F(i,j) = dt*D(i,j)*(u_old(i) - u_old(j)) +...
            dt*(B(j,i)*G(i) - B(i,j)*G(j));
    end
end

return
end