function u0 = IC_square(x)

n = length(x);
u0 = zeros(n,1);

for i = 1:n
    if ((x(i) >= 0.1) && (x(i) <= 0.3))
        u0(i) = 1;
    end
end

end