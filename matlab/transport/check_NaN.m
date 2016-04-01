function check_NaN(x)

% get size of vector
n = size(x);

% loop over elements of vector
for i = 1:n
    if (x(i) ~= x(i))
        error('NaN encountered');
    end
end

end
