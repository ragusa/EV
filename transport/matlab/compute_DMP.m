function [W_max,W_min] = compute_DMP(u_old,u_new,dt,ML,AL,b,theta,inc,periodic_BC)
% computes the upper and lower bound for max principle
%
% W_max = upper bound for max principle
% W_min = lower bound for max principle
%
% u_old = old solution
% u_new = new solution
% dt    = time step size
% ML    = lumped mass matrix
% AL    = low-order steady-state matrix; A = -(K+D)
% b     = steady-state rhs
% theta = time-discretization parameter for theta-scheme

% size of system
n = length(u_old);

W_max = zeros(n,1);
W_min = zeros(n,1);

for i = 1:n
    % compute index range of support of i
    if ~periodic_BC
        iL = max(i-1,1);
        iR = min(i+1,n);
    else
        if i == 1
            iL = n;
            iR = 2;
        elseif i == n
            iL = n-1;
            iR = 1;
        else
            iL = i-1;
            iR = i+1;
        end
    end
    % compute max and min old solution in the support of each dof
    support = [iL i iR];
    u_max_old = max(u_old(support));
    u_min_old = min(u_old(support));
    % compute max and min new solution in the support of each dof
    u_max_new = max(u_new(support));
    u_min_new = min(u_new(support));
    % compute upper and lower bounds
    aux1 = 1-(1-theta)*dt/ML(i,i)*sum(AL(i,:));
    aux2 = theta*dt/ML(i,i)*(sum(AL(i,:))-AL(i,i));
    aux3 = 1+theta*dt/ML(i,i)*AL(i,i);
    W_max(i) = (aux1*u_max_old - aux2*u_max_new + dt/ML(i,i)*b(i))/aux3;
    W_min(i) = (aux1*u_min_old - aux2*u_min_new + dt/ML(i,i)*b(i))/aux3;
end

if ~periodic_BC
    W_max(1) = inc;
    W_min(1) = inc;
end

return
end