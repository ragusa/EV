function u0 = IC_exponential_and_square(x)

n = length(x);
u0 = zeros(n,1);

for i = 1:n
    u0(i) = exp(-200*(x(i)-0.3)^2);
    if ((x(i) >= 0.6) && (x(i) <= 0.8))
        u0(i) = u0(i) + 1;
    end
end

end