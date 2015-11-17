function [Qplus,Qminus] = compute_Q(u_old,u_new,ML,W_max,W_min,AL,b,dt,theta)

n = length(u_old);

% Q vectors
Qplus = zeros(n,1);
Qminus = zeros(n,1);
for i = 1:n
    Qplus(i)  = (ML(i,i)/dt + theta*AL(i,i))*W_max(i)...
        - (ML(i,i)/dt - (1-theta)*AL(i,i))*u_old(i)...
        + AL(i,:)*(theta*u_new    + (1-theta)*u_old)...
        - AL(i,i)*(theta*u_new(i) + (1-theta)*u_old(i))...
        - b(i);
    Qminus(i) = (ML(i,i)/dt + theta*AL(i,i))*W_min(i)...
        - (ML(i,i)/dt - (1-theta)*AL(i,i))*u_old(i)...
        + AL(i,:)*(theta*u_new    + (1-theta)*u_old)...
        - AL(i,i)*(theta*u_new(i) + (1-theta)*u_old(i))...
        - b(i);
end

end