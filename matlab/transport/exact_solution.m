function u = exact_solution(x,t,initial_solution,c,sigmaL,sigmaR,qL,qR,inc,len,periodic_BC)

n = length(x);
u = zeros(1,n);

% left boundary
xL = x(1);

% middle of domain
mid = len/2;

% distance traveled
s = c*t;

u0 = @(x) inc*(x<0) + initial_solution(x)*(x>=0);

% if periodic BC, assume sigma = source = 0
if periodic_BC
    for i = 1:n
        % since IC function is not periodic, shift x(i)-s by len until
        % x(i)-s is in range (xL,xL+len), for which IC function is defined
        xs = x(i) - s;
        while xs < xL
            xs = xs + len;
        end
        u(i) = initial_solution(xs);
    end
else
    for i = 1:n
        dL = max(min(x(i),mid)-max(x(i)-s,0),0);   % distance traveled in left half
        dR = max(min(x(i),len)-max(x(i)-s,mid),0); % distance traveled in right half
        u(i) = u0(x(i)-s)*exp(-(sigmaL*dL + sigmaR*dR));
        if (sigmaL == 0)
            u(i) = u(i) + qL*dL*exp(-sigmaR*dR);
        else
            u(i) = u(i) + qL/sigmaL*(1-exp(-sigmaL*dL))*exp(-sigmaR*dR);
        end
        if (sigmaR == 0)
            u(i) = u(i) + qR*dR;
        else
            u(i) = u(i) + qR/sigmaR*(1-exp(-sigmaR*dR));
        end
    end
end

end