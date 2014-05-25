function [W_max,W_min] = compute_max_principle_bounds(u,AL,b)
% computes the upper and lower bound for max principle:
%   W_max,i = u_max,i*(1-sum(A(i,:))/A(i,i)) + b_i/A(i,i) <= u_i
%   W_min,i = u_min,i*(1-sum(A(i,:))/A(i,i)) + b_i/A(i,i) >= u_i
%
% W_max = upper bound for max principle
% W_min = lower bound for max principle
%
% u  = solution
% AL = low-order steady-state matrix; A = -(K+D)
% b  = steady-state rhs

% size of system
n = length(u);

W_max = zeros(n,1);
W_min = zeros(n,1);

for i = 1:n
    % compute index range of support of i
    i1=max(i-1,1);
    i2=min(i+1,n);
    % compute max and min old solution in the support of each dof
    u_max = max(u(i1:i2));
    u_min = min(u(i1:i2));
    % compute upper and lower bounds
    W_max(i) = u_max*(1-sum(AL(i,:))/AL(i,i)) + b(i)/AL(i,i);
    W_min(i) = u_min*(1-sum(AL(i,:))/AL(i,i)) + b(i)/AL(i,i);
end

return
end