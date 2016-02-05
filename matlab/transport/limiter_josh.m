function Flim = limiter_josh(F,Qplus,Qminus,periodic_BC)

n = size(F,1);

% compute limited fluxes
Flim = zeros(n,1);
for edge = 1:n-1
    i = edge; % upwind node
    j = i+1;  % downwind node
    if (F(i,j) > 0)
        if (i == 1)
            %Rplus_i = 1;
            Rplus_i = 0;
        else
            Rplus_i = min((Qplus(i) - Flim(i))/F(i,j),1);
        end
        Rminus_j = min((Qminus(j) - Flim(j))/F(j,i),1);
        lim = min(Rplus_i,Rminus_j);
        Flim(i) = Flim(i) + lim*F(i,j);
        Flim(j) = Flim(j) + lim*F(j,i);
    elseif (F(i,j) < 0)
        if (i == 1)
            %Rminus_i = 1;
            Rminus_i = 0;
        else
            Rminus_i = min((Qminus(i) - Flim(i))/F(i,j),1);
        end
        Rplus_j = min((Qplus(j) - Flim(j))/F(j,i),1);
        lim = min(Rminus_i,Rplus_j);
        Flim(i) = Flim(i) + lim*F(i,j);
        Flim(j) = Flim(j) + lim*F(j,i);
    end
end

end
