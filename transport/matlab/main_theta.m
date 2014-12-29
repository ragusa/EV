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
high_order_scheme = 2; % high-order scheme option
cE = 1; % coefficient for entropy residual in entropy viscosity
cJ = cE; % coefficient for jumps in entropy viscosity
entropy       = @(u) 0.5*u.^2; % entropy function
entropy_deriv = @(u) u;       % derivative of entropy function
%--------------------------------------------------------------------------
% quadrature options
%--------------------------------------------------------------------------
nq = 3; % number of quadrature points
%--------------------------------------------------------------------------
% time options
%--------------------------------------------------------------------------
theta = 0.0;     % theta for theta time integration scheme (0 = FE, 1 = BE)
CFL = 1.0;       % CFL number
ss_tol = 1.0e-5; % steady-state tolerance
n_max = 1000;    % maximum number of time steps
t_max = 10.0;    % max time to run
%--------------------------------------------------------------------------
% FCT options
%--------------------------------------------------------------------------
% DMP_option: 1 = use DMP
%             2 = use CMP
%             3 = use upwind CMP
% limiting_option: 0 = set limiting coefficients to 0 (no correction)
%                  1 = set limiting coefficients to 1 (full correction)
%                  2 = use Zalesak's limiter
%                  3 = use Josh's limiter
%
DMP_option = 1;        % DMP option
limiting_option = 2;   % limiter option
prelimit_flux = false; % option to prelimit flux
m_max = 30;            % maximum number of defect-correction iterations
dc_tol = 1e-4;         % defect-correction tolerance for discrete L2 norm
%--------------------------------------------------------------------------
% plot options
%--------------------------------------------------------------------------
compute_FCT        = false;  % compute and plot FCT solution?
compute_high_order = true;  % compute and plot high-order solution?
plot_viscosity     = false; % plot viscosities?
plot_low_order_transient  = false; % plot low-order transient?
plot_high_order_transient = true; % plot high-order transient?
plot_FCT_transient        = true; % plot FCT transient?
pausetime = 0.1;                   % time to pause for transient plots
%--------------------------------------------------------------------------
% physics options
%--------------------------------------------------------------------------
% problemID: 0: custom - use parameters below
%            1: void without source -> absorber without source
%            2: void with source    -> absorber without source
% IC_option: 0: zero
%            1: exponential pulse
%            2: exponential and square pulse
%
problemID = 1; % problem ID (if any)
IC_option = 2; % initial solution option
periodic_BC = true; % option for periodic BC; otherwise Dirichlet
len    = 1;    % length of domain
omega  = 1;    % omega
sigmaL = 0;    % sigma for left half of domain
sigmaR = 0;    % sigma for right half of domain
qL     = 0;    % source for left half of domain
qR     = 0;    % source for right half of domain
inc    = 0;    % incoming flux
speed  = 1;    % advection speed
%--------------------------------------------------------------------------
% output options
%--------------------------------------------------------------------------
save_exact_solution     = false; % option to save exact solution 
save_low_order_solution = false; % option to save low-order solution
save_FCT_solution       = false; % option to save FCT solution
%-------------------------------------------------------------------------

%% Setup

% set parameters if a problem ID was chosen
switch problemID
    case 0 % custom
        % do nothing, parameters in input section are used
    case 1 % void without source -> absorber without source
        IC_option = 0;
        periodic_BC = false;
        len    = 10;
        omega  = 1;
        sigmaL = 0;
        sigmaR = 1;
        qL     = 0;
        qR     = 0;
        inc    = 1;
        speed  = 1;
    case 2 % void with source    -> absorber without source
        IC_option = 0;
        periodic_BC = false;
        len    = 10;
        omega  = 1;
        sigmaL = 0;
        sigmaR = 1;
        qL     = 1;
        qR     = 0;
        inc    = 1;
        speed  = 1;
    case 3 % linear advection
        IC_option = 2;
        periodic_BC = true;
        len    = 1;
        omega  = 1;
        sigmaL = 0;
        sigmaR = 0;
        qL     = 0;
        qR     = 0;
        inc    = 0;
        speed  = 1;
    otherwise
        error('Invalid problem ID chosen');
end

% compute mesh quantities
n_dof = nel + 1;              % number of dofs
if periodic_BC
    n_dof = nel;
end
x = linspace(0,len,nel+1)';   % mesh points
dx_cell = diff(x);            % element sizes
dx = dx_cell(1);              % element size assuming uniform mesh
if periodic_BC
    x(end)=[];
end

% set initial conditions function and compute initial solution
switch IC_option
    case 0 % zero
        IC = @(x) zeros(length(x),1);
    case 1 % exponential pulse
        IC = @(x) exp(-200*(x-0.5).^2);
    case 2 % exponential and square pulse
        IC = @(x) exp(-200*(x-0.3).^2)+(x>=0.6).*(x<=0.8);
    otherwise
        error('Invalid IC option');
end
% compute initial solution
u0 = IC(x);

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

% compute min and max sigma and source in the support of i
sigma_upwind_min = zeros(n_dof,1);
sigma_upwind_max = zeros(n_dof,1);
q_upwind_min = zeros(n_dof,1);
q_upwind_max = zeros(n_dof,1);
for i = 1:n_dof
    i1 = max(i-1,1);
    i2 = i1;
    sigma_upwind_min(i) = min(sigma(i1:i2));
    sigma_upwind_max(i) = max(sigma(i1:i2));
    q_upwind_min(i)     = min(source(i1:i2));
    q_upwind_max(i)     = max(source(i1:i2));
end

%% Assembly

% build matrices
[MC,ML,A,AL,DL,b,viscL] = ...
    build_matrices(len,nel,omega,speed,sigma,source,...
    low_order_scheme,high_order_scheme,periodic_BC);

% check that speed is nonzero; otherwise, infinite dt will be computed
if (speed == 0)
    error('Speed cannot be zero; infinite dt will be computed');
end

% compute dt using CFL condition, even for implicit, just so that time
% steps will be equal between explicit and implicit, for comparison
dt = CFL*ML(1,1)/AL(1,1);
for i = 2:n_dof
    dt = min(dt, CFL*ML(i,i)/AL(i,i));
end

% compute transient sytem matrices and modify for Dirichlet BC
ALtr = ML + dt*theta*AL;
ALtr_mod = ALtr;
if ~periodic_BC
    ALtr_mod(1,:)=0; ALtr_mod(1,1)=1;
end
dc_preconditioner = ALtr_mod;

%% Low-order Solution

fprintf('\nComputing low-order solution...\n\n');
u_old = u0;
t = 0;
for time_step = 1:n_max
    fprintf('Time step %i: t = %f->%f',time_step,t,t+dt);
    
    % compute low-order solution
    uL = low_order_step(u_old,AL,ALtr_mod,ML,b,dt,theta,inc,periodic_BC);
    
    % test steady-state convergence
    ss_err = norm(uL-u_old,2);
    fprintf(' norm(u_new - u_old) = %e\n',ss_err);
    finished = ss_err < ss_tol || t >= t_max;
    if (finished)
        break;
    end
    
    % plot
    if (plot_low_order_transient)
        plot(x,uL);
        legend('Low-order','Location','Best');
        pause(pausetime);
    end
    
    % reset u_old
    u_old = uL;
    t = t+dt;
end
% save uL for plotting because "uL" is used in FCT loop
uL_final = uL;

%% High-order solution

if (compute_high_order)
    fprintf('\nComputing high-order solution...\n\n');
    u_old = u0;
    uH    = u0;
    t = 0;
    for time_step = 1:n_max
        fprintf('Time step %i: t = %f->%f',time_step,t,t+dt);
        
        % perform high-order step
        [uH,DH] = high_order_step(u_old,u_old,viscL,dx,x,omega,sigma,...
            source,inc,dt,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
            A,b,MC,theta,time_step,high_order_scheme,periodic_BC);
        
        % test steady-state convergence
        ss_err = norm(uH-u_old,2);
        fprintf(' norm(u_new - u_old) = %e\n',ss_err);
        finished = ss_err < ss_tol || t >= t_max;
        if (finished)
            break;
        end
        
        % plot
        if (plot_high_order_transient)
            plot(x,uH);
            legend('High-order');
            pause(pausetime);
        end
        
        % reset u_old and u_older and advance time
        u_older = u_old;
        u_old   = uH;
        t = t+dt;
    end
    
    % save uH for plotting because "uH" is used in FCT loop
    uH_final = uH;
end

%% FCT solution

if (compute_FCT)
    fprintf('\nComputing FCT solution...\n\n');
    u_old   = u0;
    u_older = u0;
    uFCT    = u0;
    t = 0;
    for time_step = 1:n_max
        fprintf('Time step %i: t = %f->%f\n',time_step,t,t+dt);
        
        % perform high-order step
        [uH,DH] = high_order_step(u_older,u_old,viscL,dx,...
            x,omega,sigma,source,inc,dt,v,dvdz,zq,wq,Jac,cE,cJ,...
            A,b,MC,theta,time_step,high_order_scheme,periodic_BC);
        
        % compute FCT solution
        %------------------------------------------------------------------
        % compute flux corrections
        F = flux_correction_matrix(u_old,uH,dt,DH,DL,MC,theta);
        
        % prelimit flux as Kuzmin suggested
        if (prelimit_flux)
            % compute low-order solution
            rhs = (ML - (1-theta)*dt*AL)*u_old + dt*b;
            if (~periodic_BC)
                rhs(1) = inc;
            end
            uL = ALtr_mod \ rhs;
            % prelimit fluxes (cancel down gradient)
            for i = 1:n_dof
                for j = max(i-1,1):min(i+1,n_dof)
                    if (F(i,j)*(u_old(i)-u_old(j)) <= 0)
                        F(i,j) = 0;
                        %fprintf('canceled for i=%i, j=%i\n',i,j);
                    end
                end
            end
        end

        % solve for new step solution using defect-correction
        dc_converged = 0; % flag for convergence of defect-correction
        for m = 1:m_max
            % compute max principle bounds
            switch DMP_option
                case 1 % DMP
                    [Wplus,Wminus] = compute_DMP(...
                        u_old,uFCT,dt,ML,AL,b,theta,inc,periodic_BC);
                case 2 % CMP
                    [Wplus,Wminus] = compute_CMP(...
                        u_old,sigma_min,sigma_max,q_min,q_max,speed*dt,inc,periodic_BC);
                case 3 % upwind CMP
                    [Wplus,Wminus] = compute_upwind_CMP(...
                        u_old,sigma_upwind_min,sigma_upwind_max,...
                        q_upwind_min,q_upwind_max,speed*dt,inc,x);
                otherwise
                    error('Invalid DMP option');
            end
            [Qplus,Qminus] = compute_Q(...
                u_old,uFCT,ML,Wplus,Wminus,AL,b,dt,theta);
            Qplus = max(Qplus,0);
            Qminus = min(Qminus,0);
            
            % compute limited fluxes
            switch limiting_option
                case 0 % no correction
                    flim = zeros(n_dof,1);
                case 1 % full correction (no limiting)
                    flim = sum(F,2);
                case 2 % Zalesak's limiter
                    flim = limiter_zalesak(F,Qplus,Qminus,periodic_BC);
                case 3 % Josh's limiter
                    flim = limiter_josh(F,Qplus,Qminus,periodic_BC);
                otherwise
                    error('Invalid limiting option');
            end
            
            % compute next FCT solution iterate using defect-correction
            rhs = (ML - (1-theta)*dt*AL)*u_old + dt*b + flim;
            if (~periodic_BC)
                rhs(1) = inc;
            end
            r = rhs - ALtr_mod*uFCT;
            du = dc_preconditioner \ r;
            uFCT = uFCT + du;
            
            % test defect-correction error
            L2norm = norm(du,2);
            fprintf('\tDC iteration %i error: %e\n',m,L2norm);
            if (L2norm < dc_tol)
                dc_converged = 1;
                break;
            end
        end
        
        % terminate if defect-correction did not converge
        if (~dc_converged)
            error('Defect correction did not converge');
        end
        
        % test steady-state convergence
        ss_err = norm(uFCT-u_old,2);
        fprintf('\tnorm(u_new - u_old) = %e\n',ss_err);
        finished = ss_err < ss_tol || t >= t_max;
        if (finished)
            break;
        end
        
        % plot
        if (plot_FCT_transient)
            plot(x,Wminus,'k:o');
            hold on;
            plot(x,Wplus,'k:x');
            plot(x,uL,'r-^');
            plot(x,uH,'b-+');
            plot(x,uFCT,'g-s');
            hold off;
            legend('W-','W+','Low','High','FCT','Location','Best');
            pause(pausetime);
        end
        
        % reset u_old and u_older and advance time
        u_older = u_old;
        u_old   = uFCT;
        t = t+dt;
    end
end

%% Plot Solution

figure(1); clf; hold on;

% plot exact solution
xx = linspace(0,len,1000)';
u_exact = exact_solution(xx,t+dt,IC,speed,sigmaL,sigmaR,qL,qR,inc,len,periodic_BC)';
plot(xx,u_exact,'k-');
legend_entries = char('Exact');

% plot low-order solution
plot(x,uL_final,'r-s');
legend_entries = char(legend_entries,'Low-order');

% plot high-order solution
if (compute_high_order)
    plot(x,uH_final,'b-+');
    legend_entries = char(legend_entries,'High-order');
end

% plot FCT solution
if (compute_FCT)
    plot(x,uFCT,'g-x');
    switch limiting_option
        case 0 % no correction
            limiter_string = 'no correction';
        case 1 % full correction (no limiting)
            limiter_string = 'not limited';
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

%% Plot viscosity

% plot viscosity if requested
if (plot_viscosity)
    figure(2); clf;
    
    % plot low-order viscosity if available
    if low_order_scheme == 2 || high_order_scheme == 2
        semilogy(x_center,viscL);
        legend_entries = char('Low-order viscosity');
    end
    
    hold on;
    
    % plot high-order viscosity if available
    if (high_order_scheme == 2)
        semilogy(x_center,viscE,'x');
        legend_entries = char(legend_entries,'Entropy viscosity');
    end
    
    % legend
    legend(legend_entries,'Location','Best');
end

%% Output

% save exact solution
if (save_exact_solution)
    exact_file = 'output/uexact.csv';
    csvwrite(exact_file,[xx,u_exact]);
end

% determine string to be appended to results for time discretization
if (theta == 0)
    time_string = 'FE';
elseif (theta == 0.5)
    time_string = 'CN';
elseif (theta == 1.0)
    time_string = 'BE';
else
    error('Invalid theta value');
end

% save low-order solution
if (save_low_order_solution)
    low_order_file = ['output/uL_',time_string,'.csv'];
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

% save FCT solution
if (save_FCT_solution)
    FCT_file = ['output/uFCT_',high_order_string,'_',time_string,'.csv'];
    csvwrite(FCT_file,[x,uFCT]);
end