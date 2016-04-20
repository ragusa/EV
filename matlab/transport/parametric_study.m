
% number of cells in coarsest mesh
n_cells_coarse = 2;

% number of refinement cycles; each subsequent cycle doubles number of cells
n_cycles = 12;

% initialize return value arrays
values1 = zeros(n_cycles,1);
values2 = zeros(n_cycles,1);

% loop over refinement cycles
n_cells = n_cells_coarse;
for i = 1:n_cycles
  % run this cycle
  [values1(i),values2(i)] = main(n_cells);

  % double number of cells
  n_cells = n_cells * 2;
end

% save convergence results
dlmwrite('output/convergence.csv',[values2,values1],',');

% Plot convergence results. For now, it is assumed that the return values are:
%   1. some quantity
%   2. mesh size h
figure;
loglog(values2, values1,'b-x');
hold on;

% plot reference slopes
slope1 = zeros(n_cycles,1);
slope2 = zeros(n_cycles,1);
slope1_end = zeros(n_cycles,1);
slope2_end = zeros(n_cycles,1);
slope1(1) = values1(1);
slope2(1) = values1(1);
slope1_end(end) = values1(end);
slope2_end(end) = values1(end);
c1 = slope1(1) / values2(1)^1;
c2 = slope2(1) / values2(1)^2;
c1_end = slope1_end(end) / values2(end)^1;
c2_end = slope2_end(end) / values2(end)^2;
for i = 2:n_cycles
  slope1(i) = c1 * values2(i)^1;
  slope2(i) = c2 * values2(i)^2;
end
for i = 1:n_cycles-1
    slope1_end(i) = c1_end * values2(i)^1;
    slope2_end(i) = c2_end * values2(i)^2;
end
loglog(values2, slope1, 'k--');
loglog(values2, slope2, 'k:');
loglog(values2, slope1_end, 'k--');
loglog(values2, slope2_end, 'k:');

% annotations
ylabel('f(h)');
xlabel('Mesh size, h');
legend('f(h)','m = 1','m = 2','Location','SouthEast');
