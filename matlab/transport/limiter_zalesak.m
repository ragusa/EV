function [flim,L] = limiter_zalesak(F,Qplus,Qminus,flim_prev,opts,...
    dirichlet_limiting_coefficient)

n = length(Qplus);

flim   = zeros(n,1);
Pplus  = zeros(n,1);
Pminus = zeros(n,1);
Rminus = zeros(n,1);
Rplus  = zeros(n,1);
L = zeros(n,n);

% P vectors
for i=1:n   
    Pplus(i) = sum(max(0,F(i,:)));
    Pminus(i)= sum(min(0,F(i,:)));
end

% R vectors
for i = 1:n
    if (Pplus(i) ~= 0)
        Rplus(i) = min(1, (Qplus(i)-flim_prev(i))/Pplus(i));
    else
        Rplus(i) = 1.0;
    end
    if (Pminus(i) ~= 0)
        Rminus(i) = min(1, (Qminus(i)-flim_prev(i))/Pminus(i));
    else
        Rminus(i) = 1.0;
    end
end

% force some bounds for limiting coefficients of dirichlet node
if (opts.impose_DirichletBC_strongly)
    Rplus(1)  = dirichlet_limiting_coefficient;
    Rminus(1) = dirichlet_limiting_coefficient;
end

% limiting coefficients
for i=1:n
    for j=1:n
        if(F(i,j)>=0)
            Lij = min(Rplus(i),Rminus(j));
            L(i,j) = Lij;
            flim(i) = flim(i) + Lij*F(i,j);
        else
            Lij = min(Rminus(i),Rplus(j));
            L(i,j) = Lij;
            flim(i) = flim(i) + Lij*F(i,j);
        end
    end
end

return
end
