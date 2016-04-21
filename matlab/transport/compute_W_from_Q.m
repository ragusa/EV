function [Wplus,Wminus] = compute_W_from_Q(u_old,u_new,ML,Qplus,Qminus,AL,b,dt,theta)

n = length(u_old);

% W vectors
Wplus  = zeros(n,1);
Wminus = zeros(n,1);
for i = 1:n
    % compute W
    Wplus(i)  = u_old(i) + (Qplus(i)  - AL(i,:)*(theta*u_new + (1-theta)*u_old) + b(i))*dt/ML(i,i);
    Wminus(i) = u_old(i) + (Qminus(i) - AL(i,:)*(theta*u_new + (1-theta)*u_old) + b(i))*dt/ML(i,i);
end

end