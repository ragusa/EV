function u = exact_solution(problemID,x,sigma,src,inc,omega,speed)

switch problemID
    case 1
        if (sigma(x) > eps)
            u = src/sigma(x)+(inc-src/sigma(x))*exp(-sigma(x)*x/omega/speed);
        else
            u = inc+src/omega/speed*x;
        end
    case 2
        n = length(x);
        u = zeros(1,n);
        for i = 1:n
            if (x(i) < 5)
                u(i) = inc;
            else
                u(i) = inc*exp(-sigma(x(i))*(x(i)-5));
            end
        end
    case 3
        u = src/sigma(x)+(inc-src/sigma(x))*exp(-sigma(x)*x/omega);
    otherwise
        error('Invalid problem ID');
end

end