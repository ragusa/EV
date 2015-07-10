function F = prelimit_fluxes(F,U)
% Prelimits the fluxes as Kuzmin suggested
%
% Output:
%  F  new flux correction matrix
%
% Input:
%  F  old flux correction matrix
%  U  solution by which to check gradient

% size of system
n = length(U);

% loop over every entry
for i = 1:n
    for j = 1:n
        if (F(i,j)*(U(i) - U(j)) <= 0)
            F(i,j) = 0;
        end
    end
end

end
