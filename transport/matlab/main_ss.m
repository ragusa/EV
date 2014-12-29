close all; clear; clc;

%% User Options
%--------------------------------------------------------------------------
% mesh options
%--------------------------------------------------------------------------
nel = 50; % number of elements
%--------------------------------------------------------------------------
% method options
%--------------------------------------------------------------------------
% low_order_scheme: 1 = algebraic low-order scheme
%                   2 = graph-theoretic low-order scheme
% high_order_scheme: 1 = Galerkin
%                    2 = Entropy viscosity
%
low_order_scheme  = 2; % low-order scheme option
high_order_scheme = 1; % high-order scheme option
cE = 1; % coefficient for entropy residual in entropy viscosity
cJ = cE; % coefficient for jumps in entropy viscosity
entropy       = @(u) 0.5*u.^2; % entropy function
entropy_deriv = @(u) u;       % derivative of entropy function
%--------------------------------------------------------------------------
% quadrature options
%--------------------------------------------------------------------------
nq = 3; % number of quadrature points
%--------------------------------------------------------------------------
% FCT options
%--------------------------------------------------------------------------
% DMP_option: 1 = use original DMP
%             2 = use new DMP
% limiting_option: 0 = set limiting coefficients to 0 (no correction)
%                  1 = set limiting coefficients to 1 (full correction)
%                  2 = use Zalesak's limiter
%                  3 = use Josh's limiter
%
compute_FCT = true;   % option to compute FCT solution
DMP_option = 2;       % DMP option
limiting_option = 2;  % limiter option
m_max = 30;           % maximum number of defect-correction iterations
dc_tol = 1e-4;        % defect-correction tolerance for discrete L2 norm
%--------------------------------------------------------------------------
% plot options
%--------------------------------------------------------------------------
plot_high_order_solution  = true; % plot high-order final solution?
%--------------------------------------------------------------------------
% physics options
%--------------------------------------------------------------------------
% problemID: 0: custom - use parameters below
%            1: void without source -> absorber without source
%            2: void with source    -> absorber without source
%
problemID = 1; % problem ID (if any)
len    = 10;   % length of domain
omega  = 1;    % omega
sigmaL = 0;    % sigma for left half of domain
sigmaR = 0;    % sigma for right half of domain
qL     = 0;    % source for left half of domain
qR     = 0;    % source for right half of domain
inc    = 1;    % incoming flux
speed  = 1;    % advection speed
%--------------------------------------------------------------------------
% output options
%--------------------------------------------------------------------------
save_exact_solution     = false; % option to save exact solution 
save_low_order_solution = false; % option to save low-order solution
%--------------------------------------------------------------------------

%% Setup

% set parameters if a problem ID was chosen
switch problemID
    case 0 % custom
        % do nothing, parameters in input section are used
    case 1 % void without source -> absorber without source
        len    = 10;
        omega  = 1;
        sigmaL = 0;
        sigmaR = 1;
        qL     = 0;
        qR     = 0;
        inc    = 1;
        speed  = 1;
    case 2 % void with source    -> absorber without source
        len    = 10;
        omega  = 1;
        sigmaL = 0;
        sigmaR = 1;
        qL     = 1;
        qR     = 0;
        inc    = 1;
        speed  = 1;
    otherwise
        error('Invalid problem ID chosen');
end

% compute mesh quantities
n_dof = nel + 1;              % number of dofs
x = linspace(0,len,nel+1)';   % mesh points
dx_cell = diff(x);            % element sizes
dx = dx_cell(1);              % element size assuming uniform mesh

% get quadrature points and weights and evaluate basis functions
[zq,wq]  = get_GL_quadrature(nq);
[v,dvdz] = get_lagrange_basis(zq,2); 
Jac = 0.5*dx_cell; % Jacobians of reference cell transformations

%% Material properties

% compute sigma and source in each cell
x_center = linspace(0.5*dx(1),len-0.5*dx(1),nel)'; % cell center positions
x_mid = len/2;                             % center of domain
sigma  = (x_center < x_mid)*sigmaL + (x_center >= x_mid)*sigmaR;
source = (x_center < x_mid)*qL     + (x_center >= x_mid)*qR;

% compute min and max sigma and source in the support of i
sigma_min = zeros(n_dof,1);
sigma_max = zeros(n_dof,1);
q_min = zeros(n_dof,1);
q_max = zeros(n_dof,1);
for i = 1:n_dof
    i1 = max(i-1,1);
    i2 = min(i,nel);
    sigma_min(i) = min(sigma(i1:i2));
    sigma_max(i) = max(sigma(i1:i2));
    q_min(i)     = min(source(i1:i2));
    q_max(i)     = max(source(i1:i2));
end

%% Assembly

% build matrices
[MC,ML,A,AL,DL,b,viscL] = ...
    build_matrices(len,nel,omega,speed,sigma,source,...
    low_order_scheme,high_order_scheme,false);

% compute sytem matrices and rhs and modify for Dirichlet BC
AL_mod = AL;
AL_mod(1,:)=0; AL_mod(1,1)=1;
b_mod = b;
b_mod(1) = inc;

%% Low-order Solution

fprintf('\tComputing low-order solution...\n');
uL = AL_mod \ b_mod;

%% High-order Solution

fprintf('\tComputing high-order solution...\n');
uH = uL;

for iter = 1:m_max
    % compute high-order diffusion matrix
    if (high_order_scheme == 2) % Entropy viscosity
        viscE = compute_entropy_viscosity_ss(...
            uH,x,omega,sigma,source,v,dvdz,wq,Jac,cE,cJ,entropy,entropy_deriv);
        DH = compute_high_order_diffusion_matrix(viscE,viscL,dx,false);
    else
        DH = spalloc(n_dof,n_dof,3*n_dof);
    end
    
    % compute high-order matrices
    AH = A + DH;
    AH_mod = AH; AH_mod(1,:)=0; AH_mod(1,1)=1;
    
    % compute residual
    res = b_mod - AH_mod*uH;
    % test convergence
    resnorm = norm(res,2);
    fprintf('\t\tIteration %i: residual error = %e\n',iter,resnorm);
    converged = resnorm < dc_tol;
    if (converged)
        fprintf('\t\tConverged at iteration %i\n',iter);
        break;
    end
    % update solution
    uH = uH + AH_mod \ res;
end

%% FCT Solution

if (compute_FCT)
    fprintf('\tComputing FCT solution...\n');
    uFCT = uL;
    
    converged = 0; % convergence flag
    for iter = 1:m_max
        % FCT solve
        %----------------
        % compute max principle bounds
        switch DMP_option
            case 1 % original DMP
                [Wplus,Wminus] = compute_DMP_ss(uFCT,AL_mod,b_mod,inc);
            case 2 % new DMP
                [Wplus,Wminus] = compute_CMP_ss(uFCT,sigma_min,sigma_max,q_min,q_max,0,inc);
            otherwise
                error('Invalid DMP option');
        end
        [Qplus,Qminus] = compute_Q_ss(uFCT,Wplus,Wminus,AL_mod,b_mod);
        
        % compute flux correction matrix
        F = flux_correction_matrix_ss(uH,DL-DH);
        
        % compute limiting coefficients
        switch limiting_option
            case 0 % Full limiter
                flim = zeros(n_dof,1);
            case 1 % No limiter
                flim = sum(F,2);
            case 2 % Zalesak limiter
                flim = limiter_zalesak(F,Qplus,Qminus,false);
            case 3 % Josh limiter
                flim = limiter_josh(F,Qplus,Qminus,false);
            otherwise
                error('Invalid limiting option');
        end
        % compute correction rhs and modify for Dirichlet BC
        rhs = b + flim;
        rhs(1) = inc;
        % compute residual
        res = rhs - AL_mod*uFCT;
        % test convergence
        resnorm = norm(res,2);
        fprintf('\t\tIteration %i: residual error = %e\n',iter,resnorm);
        converged = resnorm < dc_tol;
        if (converged)
            fprintf('\t\tConverged at iteration %i\n',iter);
            break;
        end
        % solve modified system
        uFCT = uFCT + AL_mod \ res;
    end
    
    if (~converged)
        error('\t\tSolution did not converge in %i iterations\n',m_max);
    end
end
    
%% Plot

figure(1); clf; hold on;

% plot exact solution
xx = linspace(0,len,1000);
u_exact = exact_solution_ss(xx,sigmaL,sigmaR,qL,qR,inc,len);
plot(xx,u_exact,'k');
legend_entries = char('Exact');

% plot low-order solution
plot(x,uL,'r-s');
if (exist('legend_entries','var'))
    legend_entries = char(legend_entries,'Low-order');
else
    legend_entries = char('Low-order');
end

% plot high-order solution
if (plot_high_order_solution)
    plot(x,uH,'b-+');
    legend_entries = char(legend_entries,'High-order');
end

if (compute_FCT)
    % plot FCT solution
    plot(x,uFCT,'g-x');
    switch limiting_option
        case 0 % Full limiter
            limiter_string = 'Full limiter';
        case 1 % No limiter
            limiter_string = 'No limiter';
        case 2 % Zalesak limiter
            limiter_string = 'Zalesak limiter';
        case 3 % Josh limiter
            limiter_string = 'Josh limiter';
        otherwise
            error('Invalid limiting option');
    end
    FCT_legend_string = ['FCT, ',limiter_string];
    legend_entries = char(legend_entries,FCT_legend_string);
end

% legend
legend(legend_entries,'Location','Best');