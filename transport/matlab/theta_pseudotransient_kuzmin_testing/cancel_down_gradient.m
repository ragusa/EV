function F = cancel_down_gradient(F,u_aux)

n = length(u_aux);

for i = 1:n
    for j = 1:n
        if (F(i,j)*(u_aux(i) - u_aux(j)) < 0)
            F(i,j) = 0;
        end
    end
end

end