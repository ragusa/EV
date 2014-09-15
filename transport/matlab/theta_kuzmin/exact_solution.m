function u = exact_solution(IC,x,t,len)
% assumes domain is of length 1

n = length(x);

u = zeros(n,1);
for i = 1:n
    xx = x(i) - t;
    while (xx < 0)
        xx = xx + len;
    end
    u(i) = IC(xx);
end

end