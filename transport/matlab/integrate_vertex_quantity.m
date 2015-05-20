% Computes the integral of some vertex quantity over the domain;
% assumes quantity to be piecewise linear.
function integral = integrate_vertex_quantity(x,y)

% compute quantity at cell centers
y_cell = 0.5*(y(1:end-1) + y(2:end));

% compute cell widths
dx = diff(x);

% compute integral, which is simply a sum
integral = sum(y_cell.*dx);

end