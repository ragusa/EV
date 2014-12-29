function [b,dbdz] = get_lagrange_basis(uq,dofs_per_cell)
% computes Lagrangian basis functions and their first derivatives
%  for any order p

% number of quadrature points
nq = length(uq);

% Lagrangian basis functions are equispaced between (-1,1)
xd = linspace(-1,1,dofs_per_cell);

% compute basis function values at each quadrature point
b = zeros(dofs_per_cell,nq);
for i = 1:dofs_per_cell
    numerator = 1.0;
    denominator = 1.0;
    for j = 1:dofs_per_cell
        if (j ~= i)
            numerator = numerator .* (uq'   - xd(j));
            denominator = denominator .* (xd(i) - xd(j));
        end
    end
    b(i,:) = numerator ./ denominator;
end

% compute basis function derivatives at each quadrature point
dbdz = zeros(dofs_per_cell,nq);
for i = 1:dofs_per_cell
    sum = 0.0;
    denominator = 1.0;
    for j = 1:dofs_per_cell
        if (j ~= i)
            numerator = 1.0;
            for k = 1:dofs_per_cell
                if ((k ~= i)&&(k ~= j))
                    numerator = numerator .* (uq' - xd(k));
                end
            end
            sum = sum + numerator;
            denominator = denominator .* (xd(i) - xd(j));
        end
    end
    dbdz(i,:) = sum ./ denominator;
end

return
end
