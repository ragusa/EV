close all; clear; clc;

%% User Options
%--------------------------------------------------------------------------
% mesh options
%--------------------------------------------------------------------------
nel = 10; % number of elements
%--------------------------------------------------------------------------
% method options
%--------------------------------------------------------------------------
% low_order_scheme: 1 = algebraic low-order scheme
%                   2 = graph-theoretic low-order scheme
% high_order_scheme: 1 = Galerkin
%                    2 = Entropy viscosity
%
low_order_scheme  = 2; % low-order scheme option
high_order_scheme = 2; % high-order scheme option
cE = 1.0; % coefficient for entropy residual in entropy viscosity
cJ = 1.0; % coefficient for jumps in entropy viscosity
entropy       = @(u) 0.5*u.^2; % entropy function
entropy_deriv = @(u) u;       % derivative of entropy function
impose_DirichletBC_strongly = true; % impose Dirichlet BC strongly?
%--------------------------------------------------------------------------
% quadrature options
%--------------------------------------------------------------------------
nq = 3; % number of quadrature points
%--------------------------------------------------------------------------
% FCT options
%--------------------------------------------------------------------------
% DMP_option: 1 = DMP
%             2 = max/min(DMP,CMP)
% limiting_option: 0 = set limiting coefficients to 0 (no correction)
%                  1 = set limiting coefficients to 1 (full correction)
%                  2 = use Zalesak's limiter
%                  3 = use Josh's limiter
%
compute_FCT = false;   % option to compute FCT solution
DMP_option = 2;       % DMP option
limiting_option = 2;  % limiter option
m_max = 100;           % maximum number of defect-correction iterations
dc_tol = 1e-10;        % defect-correction tolerance for discrete L2 norm
%--------------------------------------------------------------------------
% physics options
%--------------------------------------------------------------------------
% problemID: 0: custom - use parameters below
%            1: void without source -> absorber without source
%            2: void with source    -> absorber without source
%
problemID = 3; % problem ID (if any)
x_min = 0.0;   % left end of domain
x_max = 1.0;   % right end of domain
mu     = 1;    % cos(angle)
sigma  = @(x) 1;   % cross section function
source = @(x,t) 0; % source function
inc    = 1;    % incoming flux
speed  = 1;    % advection speed
%--------------------------------------------------------------------------
% output options
%--------------------------------------------------------------------------
plot_FCT_iterate = true; % option to plot FCT and bounds during iteration
pausetime = 0.0; % time to pause if plotting during iteration
plot_bounds = true; % option to plot maximum principle bounds
save_exact_solution      = false; % option to save exact solution 
save_low_order_solution  = false; % option to save low-order solution
save_high_order_solution = true; % option to save high-order solution
save_FCT_solution        = false; % option to save FCT solution
%--------------------------------------------------------------------------

%% Setup

% set parameters if a problem ID was chosen
switch problemID
    case 0 % custom
        % do nothing, parameters in input section are used
    case 1 % pure absorber without source
        x_min = 0.0;
        x_max = 1.0;
        inc    = 1.0;
        mu     = 1.0;
        sigma  = @(x) 10.0;
        source = @(x,t) 0.0;
        speed  = 1;
        exact_solution_known = true;
        exact = @(x) exp(-10*x);
    case 2 % void without source -> absorber without source
        x_min = 0.0;
        x_max = 1.0;
        inc    = 1.0;
        mu     = 1.0;
        sigma  = @(x) 10.0*(x >= 0.5);
        source = @(x,t) 0.0;
        speed  = 1;
        exact_solution_known = true;
        exact = @(x) (x<0.5) + (x>=0.5).*(exp(-10*(x-0.5)));
    case 3 % void with source -> absorber without source
        x_min = 0.0;
        x_max = 1.0;
        inc    = 0.0;
        mu     = 1.0;
        sigma  = @(x) 10.0*(x >= 0.5);
        source = @(x,t) 1.0*(x < 0.5);
        speed  = 1;
        exact_solution_known = true;
        exact = @(x) x.*(x<0.5) + 0.5*exp(-10*(x-0.5)).*(x>=0.5);
    case 4
        x_min = 0.0;
        x_max = 1.0;
        inc    = 0.0;
        mu     = 1.0;
        sigma  = @(x) 10.0;
        source = @(x,t) 1.0*(x < 0.5);
        speed  = 1;
        exact_solution_known = false;
        exact = @(x) x.*(x<0.5) + 0.5*exp(-10*(x-0.5)).*(x>=0.5);
    case 5 % MMS-1
        x_min = 0.0;
        x_max = 1.0;
        inc    = 1.0;
        mu     = 1.0;
        sigma  = @(x) 1.0;
        source = @(x,t) pi*cos(pi*x)+sin(pi*x)+1;
        speed  = 1;
        exact_solution_known = true;
        exact = @(x,t) sin(pi*x)+1;
    otherwise
        error('Invalid problem ID chosen');
end

% compute mesh quantities
n_dof = nel + 1;                  % number of dofs
x = linspace(x_min,x_max,nel+1)'; % mesh points
dx = diff(x);                     % element sizes

% get quadrature points and weights and evaluate basis functions
[zq,wq]  = get_GL_quadrature(nq);
[v,dvdz] = get_lagrange_basis(zq,2); 
Jac = 0.5*dx; % Jacobians of reference cell transformations

%% Material properties

% compute min and max sigma and source in the support of i
sigma_min  = 1e15*ones(n_dof,1);
sigma_max  = zeros(n_dof,1);
source_min = 1e15*ones(n_dof,1);
source_max = zeros(n_dof,1);

%% Assembly

% create connectivity array
connectivity = [linspace(1,nel,nel)' linspace(2,nel+1,nel)'];

% assemble mass matrices, inviscid and low-order steady-state matrices,
% low-order viscosity and corresponding artificial diffusion matrix
[MC,ML,A,AL,DL,viscL] = build_matrices(...
    nq,zq,wq,v,dvdz,Jac,x,dx,nel,n_dof,connectivity,...
    mu,speed,sigma,low_order_scheme,~impose_DirichletBC_strongly);

% assemble steady-state rhs at time 0
b = assemble_ss_rhs(nq,zq,wq,v,Jac,x,nel,n_dof,connectivity,...
    source,mu,speed,0,inc,~impose_DirichletBC_strongly);

% compute sytem matrices and rhs and modify for Dirichlet BC
AL_mod = AL;
b_mod = b;
if (impose_DirichletBC_strongly)
    AL_mod(1,:)=0; AL_mod(1,1)=1;
    b_mod(1) = inc;
end

%% Low-order Solution

fprintf('\tComputing low-order solution...\n');
uL = AL_mod \ b_mod;

%% High-order Solution

fprintf('\tComputing high-order solution...\n');
uH = uL;

for iter = 1:m_max
    % compute high-order diffusion matrix
    if (high_order_scheme == 2) % Entropy viscosity
        viscE = compute_entropy_viscosity(...
            uH,uH,x,mu,sigma,source,1.0,...
            v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,false);
        DH = compute_high_order_diffusion_matrix(viscE,viscL,dx,false);
    else
        DH = spalloc(n_dof,n_dof,3*n_dof);
    end
    
    % compute high-order matrices
    AH = A + DH;
    AH_mod = AH;
    if (impose_DirichletBC_strongly)
        AH_mod(1,:)=0; AH_mod(1,1)=1;
    end
    
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

% report if the solution did not converge
if (~converged)
    error('\t\tHigh-order Solution did not converge in %i iterations\n',m_max);
end

%% FCT Solution

if (compute_FCT)
    fprintf('\tComputing FCT solution...\n');
            
    % compute flux correction matrix
    F = flux_correction_matrix_ss(uH,DL-DH);
    
    uFCT = uL;
    converged = 0; % convergence flag
    for iter = 1:m_max
        % FCT solve
        %----------------
        % compute max principle bounds
        [Wplus,Wminus] = compute_DMP_ss(uFCT,AL_mod,b_mod,inc);
        if DMP_option == 2
            [WplusCMP,WminusCMP] = compute_CMP_ss(uFCT,sigma_min,...
                    sigma_max,source_min,source_max,0,inc);
            Wplus = max(Wplus,WplusCMP);
            Wminus = min(Wminus,WminusCMP);
        end
        [Qplus,Qminus] = compute_Q_ss(uFCT,Wplus,Wminus,AL_mod,b_mod);
        
        % plot iterate of FCT solution
        if (plot_FCT_iterate)
            clf;
            hold on;
            plot(x,uH,'b-+');
            plot(x,uFCT,'g-x');
            plot(x,Wminus,'k--');
            plot(x,Wplus,'k--');
            hold off;
            legend('High-order',['FCT,l=',num2str(iter)],['W-,l=',num2str(iter)],['W+,l=',num2str(iter)]);
            pause(pausetime);
        end
        
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
        if (impose_DirichletBC_strongly)
            rhs(1) = inc;
        end
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
        error('\t\tFCT Solution did not converge in %i iterations\n',m_max);
    end
end
    
%% Plot

figure(1); clf; hold on;

% plot exact solution
if (exact_solution_known)
    xx = linspace(x_min,x_max,1000);
    u_exact = exact(xx);
    plot(xx,u_exact,'k');
    legend_entries = char('Exact');
end

% plot low-order solution
plot(x,uL,'r-s');
if (exact_solution_known)
    legend_entries = char(legend_entries,'Low-order');
else
    legend_entries = char('Low-order');
end

% plot high-order solution
plot(x,uH,'b-+');
legend_entries = char(legend_entries,'High-order');

if (compute_FCT)
    % plot FCT solution
    plot(x,uFCT,'g-x');
    legend_entries = char(legend_entries,'FCT');
    
    % plot bounds
    if (plot_bounds)
        plot(x,Wminus,'k--');
        plot(x,Wplus,'k--');
        legend_entries = char(legend_entries,'W-','W+');
    end
end

% legend
legend(legend_entries,'Location','Best');

%% Output

% save exact solution
if (save_exact_solution)
    exact_file = 'output/uexact.csv';
    csvwrite(exact_file,[xx,u_exact]);
end

% save low-order solution
if (save_low_order_solution)
    low_order_file = ['output/uL_ss.csv'];
    csvwrite(low_order_file,[x,uL]);
end

% determine string to be appended for high-order scheme
switch high_order_scheme
    case 1
        high_order_string = 'Gal';
    case 2
        high_order_string = 'EV';
    otherwise
        error('Invalid high-order scheme chosen');
end

% save high-order solution
if (save_high_order_solution)
    high_order_file = ['output/uH_',high_order_string,'_ss.csv'];
    csvwrite(high_order_file,[x,uH]);
end

% save FCT solution
if (save_FCT_solution)
    if (compute_FCT)
        FCT_file = ['output/uFCT_',high_order_string,'_ss.csv'];
        csvwrite(FCT_file,[x,uFCT]);
    end
end