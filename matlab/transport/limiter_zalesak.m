function flim = limiter_zalesak(F,Qplus,Qminus,periodic_BC,...
    dirichlet_limiting_coefficient)

n = length(Qplus);

flim   = zeros(n,1);
Pplus  = zeros(n,1);
Pminus = zeros(n,1);
Rminus = zeros(n,1);
Rplus  = zeros(n,1);

% P vectors
for i=1:n   
    Pplus(i) = sum(max(0,F(i,:)));
    Pminus(i)= sum(min(0,F(i,:)));
end

% R vectors
for i = 1:n
    if (Pplus(i) ~= 0)
        Rplus(i) = min(1, Qplus(i)/Pplus(i));
    else
        Rplus(i) = 1.0;
    end
    if (Pminus(i) ~= 0)
        Rminus(i) = min(1, Qminus(i)/Pminus(i));
    else
        Rminus(i) = 1.0;
    end
end

% force some bounds for limiting coefficients of dirichlet node
if ~periodic_BC
    Rplus(1)  = dirichlet_limiting_coefficient;
    Rminus(1) = dirichlet_limiting_coefficient;
end

% limiting coefficients
for i=1:n
    for j=1:n
        if(F(i,j)>=0)
            Lij = min(Rplus(i),Rminus(j));
            flim(i) = flim(i) + Lij*F(i,j);
        else
            Lij = min(Rminus(i),Rplus(j));
            flim(i) = flim(i) + Lij*F(i,j);
        end
    end
end

return
end
