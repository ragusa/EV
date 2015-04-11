function F = flux_correction_matrix_ss(u,D)
% computes Kuzmin's flux correction matrix
%
% F     = flux correction matrix
%
% uH    = high-order solution
% D     = Kuzmin's artificial dissipation matrix

% size of system
n = length(u);

% compute flux correction matrix
F = zeros(n,n);

for i = 1:n
    for j = 1:n
        F(i,j) = D(i,j)*(u(j)-u(i));
    end
end

return
end