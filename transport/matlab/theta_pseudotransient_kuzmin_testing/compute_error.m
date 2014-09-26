function err = compute_error(u,xpoints)

% number of elements
nel = length(xpoints) - 1;

err = 0;
% loop over elements
for iel=1:nel
    x1 = xpoints(iel); x2 = xpoints(iel+1);
    tol = 1e-6;
    err = err + quad(@(x)my_integrand(x,x1,x2,u(iel:iel+1)),x1,x2,tol);    
end
 
err = sqrt(err);

end

%%

function out = my_integrand(x,x1,x2,u)

out = (local_solution(x,x1,x2,u)-exp(-x)).^2;

end

%%

function uloc_x = local_solution(x,x1,x2,uloc)

uloc_x = uloc(1)*b1(x,x1,x2) + uloc(2)*b2(x,x1,x2);

end

%%

function b = b1(x,x1,x2)
b = (x-x2)/(x1-x2);
end

%%

function b = b2(x,x1,x2)
b = (x1-x)/(x1-x2);
end

