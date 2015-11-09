function flim = limiter_zalesak(F,Qplus,Qminus,periodic_BC)

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

if ~ periodic_BC
    % set R+ and R- = 1 for Dirichlet node as Kuzmin recommended. This
    % prevents L_ij from automatically being 0 for j in the support of i
    Rplus(1)  = 1.0;
    Rminus(1) = 1.0;
end

% limiting coefficients
for i=1:n
    for j=1:n
        if(F(i,j)>=0)
            flim(i) = flim(i) + min(Rplus(i),Rminus(j))*F(i,j);
        else
            flim(i) = flim(i) + min(Rminus(i),Rplus(j))*F(i,j);
        end
    end
end

return
end