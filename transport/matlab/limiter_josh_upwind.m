function Flim = limiter_josh_upwind(F,Qplus,Qminus)

n = size(F,1);

% compute limited fluxes
Flim = zeros(n,1);
for edge = 1:n-1
    i = edge; % upwind node
    j = i+1;  % downwind node
    flux = F(i,j);
    if (flux > 0)
        if (Flim(i) > Qplus(i))
            flux = Qplus(i) - Flim(i);
            lim = 1;
        elseif (Flim(i) < Qminus(i))
            lim_min = (Qminus(i) - Flim(i))/flux;
            lim_max = (Qplus(i)  - Flim(i))/flux;
            lim = lim_min + max(0,min(lim_max,1)-lim_min);
        else
            lim = min((Qplus(i)-Flim(i))/flux, 1);
        end
    elseif (flux < 0)
        if (Flim(i) < Qminus(i))
            flux = Qminus(i) - Flim(i);
            lim = 1;
        elseif (Flim(i) > Qplus(i))
            lim_min = (Qplus(i)  - Flim(i))/flux;
            lim_max = (Qminus(i) - Flim(i))/flux;
            lim = lim_min + max(0,min(lim_max,1)-lim_min);
        else
            lim = min((Qminus(i)-Flim(i))/flux, 1);
        end
    else
        lim = 0;
    end
    Flim(i) = Flim(i) + lim*flux;
    Flim(j) = Flim(j) + lim*flux;
end

end