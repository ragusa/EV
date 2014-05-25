function [W_max,W_min] = compute_max_principle_bounds(u_old,u_new,dt,ML,AL,b,use_explicit)
% computes the upper and lower bound for max principle:
%   W_max,i = u_old,min,i*(1-dt/m_i*sum(A(i,:)) + dt/m_i*b_i <= u_new,i
%   W_min,i = u_old,max,i*(1-dt/m_i*sum(A(i,:)) + dt/m_i*b_i >= u_new,i
%
% W_max = upper bound for max principle
% W_min = lower bound for max principle
%
% u_old = old solution
% dt    = time step size
% ML    = lumped mass matrix
% AL    = low-order steady-state matrix; A = -(K+D)
% b     = steady-state rhs

% size of system
n = length(u_old);

W_max = zeros(n,1);
W_min = zeros(n,1);

if (use_explicit) % explicit Euler
    for i = 1:n
        % compute index range of support of i
        i1=max(i-1,1);
        i2=min(i+1,n);
        % compute max and min old solution in the support of each dof
        u_max = max(u_old(i1:i2));
        u_min = min(u_old(i1:i2));
        % compute upper and lower bounds
        W_max(i) = u_max*(1-dt/ML(i,i)*sum(AL(i,:))) + dt/ML(i,i)*b(i);
        W_min(i) = u_min*(1-dt/ML(i,i)*sum(AL(i,:))) + dt/ML(i,i)*b(i);
    end
else              % implicit Euler
    for i = 1:n
        % compute index range of support of i
        i1=max(i-1,1);
        i2=min(i+1,n);
        % compute max and min old solution in the support of each dof
        u_max = max(u_new(i1:i2));
        u_min = min(u_new(i1:i2));
        % compute upper and lower bounds
        aux = dt/ML(i,i)*(sum(AL(i,:))-AL(i,i));
        W_max(i) = (u_old(i) - aux*u_max + dt/ML(i,i)*b(i))/(1+dt/ML(i,i)*AL(i,i));
        W_min(i) = (u_old(i) - aux*u_min + dt/ML(i,i)*b(i))/(1+dt/ML(i,i)*AL(i,i));
    end
end

return
end