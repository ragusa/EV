
% read data
uL_data = csvread('output/uL_ss.csv');
P_matrix = csvread('output/P.csv');
AL = csvread('output/AL.csv');
uexact_data = csvread('output/uexact.csv');

% strongly impose Dirichlet BC?
strong_dirichlet = true;

% dirichlet value
inc = 0.0;

% extract columns from data
x = uL_data(:,1);
uL = uL_data(:,2);
xexact = uexact_data(:,1);
uexact = uexact_data(:,2);

% convert flux matrix to flux vector
P = convert_matrix_to_edge_vector(P_matrix);

% call GUI
manual_limiting_ss_gui(x,uL,AL,inc,P,xexact,uexact,strong_dirichlet);