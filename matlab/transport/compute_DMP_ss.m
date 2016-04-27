function [Wplus,Wminus] = compute_DMP_ss(u,AL,b,inc,opts)

n = length(u);

Wplus = zeros(n,1);
Wminus = zeros(n,1);

for i = 1:n
    % compute index range of support of i
    i1 = max(i-1,1);
    i2 = min(i+1,n);
    % compute max and min old solution in the support of each dof
    u_max = -1.0e15;
    u_min = 1.0e15;
    for j = i1:i2
      if (j ~= i)
        u_max = max(u_max,u(j));
        u_min = min(u_min,u(j));
      end
    end
       
    % compute bounds
    Wplus(i)  = -(sum(AL(i,:))-AL(i,i))/AL(i,i)*u_max + b(i)/AL(i,i);
    Wminus(i) = -(sum(AL(i,:))-AL(i,i))/AL(i,i)*u_min + b(i)/AL(i,i);
end

% if using strong Dirichlet BC, then no DMP applies
if (opts.modify_for_strong_DirichletBC)
    Wplus(1)  = inc;
    Wminus(1) = inc;
else
    if (AL(1,1) < 0.0)
        error('DMP does not apply because AL(1,1) < 0');
    end
end

end
