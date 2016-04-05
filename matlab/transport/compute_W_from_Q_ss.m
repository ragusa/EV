function [Wplus,Wminus] = compute_W_from_Q_ss(u,Qplus,Qminus,AL,b)

n = length(u);

% W vectors
Wplus  = zeros(n,1);
Wminus = zeros(n,1);
for i = 1:n
    % compute W
    Wplus(i)   = (Qplus(i)  + b(i) - AL(i,:)*u + AL(i,i)*u(i)) / AL(i,i);
    Wminus(i)  = (Qminus(i) + b(i) - AL(i,:)*u + AL(i,i)*u(i)) / AL(i,i);
end

end