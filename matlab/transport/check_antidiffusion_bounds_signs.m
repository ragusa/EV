function check_antidiffusion_bounds_signs(Qplus,Qminus)

% size of vector
n = length(Qplus);

small = 1.0e-12;

% skip the first node because it's assumed to be Dirichlet
for i = 2:n
  % check upper bound
  if (Qplus(i) < -small)
    error('Q+(%i) < 0: Q+(%i) = %f',i,i,Qplus(i));
  end
  % check lower bound
  if (Qminus(i) > small)
    error('Q-(%i) > 0: Q-(%i) = %f',i,i,Qminus(i));
  end
end

end
