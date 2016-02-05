function vec = convert_matrix_to_edge_vector(A)

n = size(A)(1) - 1;

vec = zeros(n,1);

for i = 1:n
    vec(i) = A(i+1,i);
end

end
