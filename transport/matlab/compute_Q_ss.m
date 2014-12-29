function [Qplus,Qminus] = compute_Q_ss(u,Wplus,Wminus,AL,b)

n = length(u);

% Q vectors
Qplus  = zeros(n,1);
Qminus = zeros(n,1);
for i = 1:n
    % compute Q
    Qplus(i)  = AL(i,i)*Wplus(i)  + AL(i,:)*u - AL(i,i)*u(i) - b(i);
    Qminus(i) = AL(i,i)*Wminus(i) + AL(i,:)*u - AL(i,i)*u(i) - b(i);
end

end