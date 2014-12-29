function u = exact_solution_ss(x,sigmaL,sigmaR,qL,qR,inc,len)
% Solves the PDE
%
% omega*dphi/dx + sigma(x)*psi(x) = q(x)
% phi(0) = inc

n = length(x);
u = zeros(1,n);

% middle of domain
mid = len/2;

for i = 1:n
    if (x(i) <= mid) % left half of domain
        if (sigmaL ~= 0)
            u(i) = inc*exp(-sigmaL*x(i)) + qL/sigmaL*(1-exp(-sigmaL*x(i)));
        else
            u(i) = inc*exp(-sigmaL*x(i)) + qL*x(i);
        end
    else % right half of domain
        % compute solution at interface
        if (sigmaL ~= 0)
            u_mid = inc*exp(-sigmaL*mid) + qL/sigmaL*(1-exp(-sigmaL*mid));
        else
            u_mid = inc*exp(-sigmaL*mid) + qL*mid;
        end
        % compute solution at node
        if (sigmaR ~= 0)
            u(i) = u_mid*exp(-sigmaR*(x(i)-mid)) + qR/sigmaR*(1-exp(-sigmaR*(x(i)-mid)));
        else
            u(i) = u_mid*exp(-sigmaR*(x(i)-mid)) + qR*(x(i)-mid);
        end
    end
end

end