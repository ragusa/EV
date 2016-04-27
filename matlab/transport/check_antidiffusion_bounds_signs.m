function check_antidiffusion_bounds_signs(Qplus,Qminus)

% size of vector
n = length(Qplus);

small = 1.0e-12;

% may need to skip first node when using strong Dirichlet BC
for i = 1:n
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
