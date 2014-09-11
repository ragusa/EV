function Flim = enforce_LED(Flim,u_aux)

n = length(u_aux);

for i = 1:n
    % compute index range of support of i
    i1=max(i-1,1);
    i2=min(i+1,n);
    % compute max and min old solution in the support of each dof
    u_max = max(u_aux(i1:i2));
    u_min = min(u_aux(i1:i2));
    
    if (u_aux(i) == u_max)
        for j = 1:n
            if (Flim(i,j) > 0)
                Flim(i,j) = 0;
            end
        end
    end
    if (u_aux(i) == u_min)
        for j = 1:n
            if (Flim(i,j) < 0)
                Flim(i,j) = 0;
            end
        end
    end
end


end