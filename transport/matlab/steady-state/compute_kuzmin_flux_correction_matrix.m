function F = compute_kuzmin_flux_correction_matrix(uH,D)
% computes Kuzmin's flux correction matrix
%
% F     = flux correction matrix
%
% uH    = high-order solution
% D     = Kuzmin's artificial dissipation matrix

% size of system
n = length(uH);

% compute flux correction matrix
F = zeros(n,n);
for i = 1:n
    for j = 1:n
        F(i,j) = D(i,j)*(uH(i)-uH(j));
    end
end

return
end