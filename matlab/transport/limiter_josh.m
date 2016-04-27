function Flim = limiter_josh(F,Qplus,Qminus,opts,...
    dirichlet_limiting_coefficient)            

n = size(F,1);

% compute limited fluxes
Flim = zeros(n,1);
for edge = 1:n-1
    i = edge; % upwind node
    j = i+1;  % downwind node
    if (F(i,j) > 0)
        Rplus_i = min((Qplus(i) - Flim(i))/F(i,j),1);
        if (opts.impose_DirichletBC_strongly && i == 1)
            Rplus_i = dirichlet_limiting_coefficient;
        end
        Rminus_j = min((Qminus(j) - Flim(j))/F(j,i),1);
        lim = min(Rplus_i,Rminus_j);
        Flim(i) = Flim(i) + lim*F(i,j);
        Flim(j) = Flim(j) + lim*F(j,i);
    elseif (F(i,j) < 0)
        Rminus_i = min((Qminus(i) - Flim(i))/F(i,j),1);
        if (opts.impose_DirichletBC_strongly && i == 1)
            Rminus_i = dirichlet_limiting_coefficient;
        end
        Rplus_j = min((Qplus(j) - Flim(j))/F(j,i),1);
        lim = min(Rminus_i,Rplus_j);
        Flim(i) = Flim(i) + lim*F(i,j);
        Flim(j) = Flim(j) + lim*F(j,i);
    end
end

end
