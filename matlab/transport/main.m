close all; clear; clc;

%% User Options
%--------------------------------------------------------------------------
% finite element options
%--------------------------------------------------------------------------
mesh.n_cell = 2^5;                   % number of elements
impose_DirichletBC_strongly = true; % impose Dirichlet BC strongly?
quadrature.nq = 3;                  % number of quadrature points per cell
%--------------------------------------------------------------------------
% spatial method options
%--------------------------------------------------------------------------
compute_low_order  = true; % compute and plot low-order solution?
compute_high_order = true; % compute and plot high-order solution?
compute_FCT        = true; % compute and plot FCT solution?

% low_order_scheme: 1 = algebraic low-order scheme
%                   2 = graph-theoretic low-order scheme
low_order_scheme  = 2;

% high_order_scheme: 1 = Galerkin
%                    2 = Entropy viscosity
high_order_scheme = 1;

% entropy viscosity options:
ev.cE = 0.1; % coefficient for entropy residual in entropy viscosity
ev.cJ = ev.cE*1; % coefficient for jumps in entropy viscosity
ev.entropy       = @(u) 0.5*u.^2; % entropy function
ev.entropy_deriv = @(u) u;        % derivative of entropy function
ev.smooth_entropy_viscosity = false; % option to smooth entropy viscosity
ev.smoothing_weight = 0.0; % weight for center value in smoothing
%--------------------------------------------------------------------------
% time options
%--------------------------------------------------------------------------
% temporal scheme: 0 = steady-state
%                  1 = SSPRK(1,1) (Explicit Euler)
%                  2 = SSPRK(3,3) (Shu-Osher)
%                  3 = theta method
temporal_scheme = 0; % temporal discretization scheme

theta = 0.5;     % theta parameter to use if using a theta method
CFL = 0.5;       % CFL number
ss_tol = 1.0e-5; % steady-state tolerance
t_end = 0.3;     % max time to run
%--------------------------------------------------------------------------
% FCT options
%--------------------------------------------------------------------------
% DMP option: 1 = low-order DMP
%             2 = max/min(low-order DMP, analytic DMP)
DMP_option = 2;

% limiter option: 0 = All 0 (no correction; low-order)
%                 1 = All 1 (full correction; high-order)
%                 2 = Zalesak limiter
%                 3 = Josh limiter
limiting_option = 3;

% FCT initialization option: 1 = zeros
%                            2 = low-order solution
%                            3 = high-order solution
FCT_initialization = 3;

% prelimit correction fluxes: 0 = do not prelimit, 1 = prelimit
prelimit = 0;
%--------------------------------------------------------------------------
% physics options
%--------------------------------------------------------------------------
% problemID: 0: custom - use parameters below
%            1: pure absorber without source
%            2: void without source -> absorber without source
%            3: void with    source -> absorber without source
%            4: void
%            5: MMS-1
%            6: MMS: TR: u = x*t, SS: u = x
problemID = 6;

% IC_option: 0: zero
%            1: exponential pulse
%            2: exponential and square pulse
IC_option = 0;

impose_BC_on_IC = 1; % option to impose Dirichlet BC on IC

mesh.x_min = 0.0;         % left end of domain
mesh.x_max = 1.0;         % right end of domain
phys.periodic_BC = false; % option for periodic BC; otherwise Dirichlet
phys.mu     = 1;          % cos(angle)
phys.sigma  = @(x,t) 0;   % cross section function
phys.source = @(x,t) 0;   % source function
phys.inc    = 1;          % incoming flux
phys.speed  = 1;          % advection speed
source_is_time_dependent = false; % is source time-dependent?
%--------------------------------------------------------------------------
% nonlinear solver options
%--------------------------------------------------------------------------
max_iter = 100;           % maximum number of nonlinear solver iterations
nonlin_tol = 1e-10;       % nonlinear solver tolerance for discrete L2 norm
relaxation_parameter = 1.0; % relaxation parameter for iteration
%--------------------------------------------------------------------------
% plot options
%--------------------------------------------------------------------------
plot_viscosity            = false; % plot viscosities?
plot_low_order_transient  = false; % plot low-order transient?
plot_high_order_transient = false; % plot high-order transient?
plot_FCT_transient        = false; % plot FCT transient?
pausetime = 0.01;                   % time to pause for transient plots
legend_location           = 'NorthEast'; % location of plot legend
%--------------------------------------------------------------------------
% output options
%--------------------------------------------------------------------------
save_exact_solution      = false; % option to save exact solution 
save_low_order_solution  = false; % option to save low-order solution
save_high_order_solution = false; % option to save high-order solution
save_FCT_solution        = false; % option to save FCT solution
%-------------------------------------------------------------------------

%% Define Problem

% set parameters if a problem ID was chosen
switch problemID
    case 0 % custom
        % do nothing, parameters in input section are used
    case 1 % pure absorber without source
        mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 1.0;
        phys.mu     = 1.0;
        phys.sigma  = @(x,t) 1.0;
        phys.source = @(x,t) 0.0;
        phys.speed  = 1;
        
        IC_option = 0;
        source_is_time_dependent = false;
        exact_solution_known = true;
        if temporal_scheme == 0
          exact = @(x,t) exp(-1*x);
        else
          exact = @(x,t) (t>=x).*exp(-10*x);
        end
    case 2 % void without source -> absorber without source
        mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 1.0;
        phys.mu     = 1.0;
        phys.sigma  = @(x,t) 10.0*(x >= 0.5);
        phys.source = @(x,t) 0.0;
        phys.speed  = 1;
        
        IC_option = 0;
        source_is_time_dependent = false;
        exact_solution_known = true;
        if temporal_scheme == 0
          exact = @(x,t) (x<0.5) + (x>=0.5).*(exp(-10*(x-0.5)));
        else
          exact = @(x,t) (t>=x).*((x<0.5) + (x>=0.5).*(exp(-10*(x-0.5))));
        end
    case 3 % void with source -> absorber without source
        mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 0.0;
        phys.mu     = 1.0;
        sigma_value = 10.0;
        phys.sigma  = @(x,t) sigma_value*(x >= 0.5);
        source_value = 1.0;
        phys.source = @(x,t) source_value*(x < 0.5);
        phys.speed  = 1;
        
        IC_option = 0;
        source_is_time_dependent = false;
        exact_solution_known = true;
        if temporal_scheme == 0
            exact = @(x,t) x.*(x<0.5) + 0.5*exp(-10*(x-0.5)).*(x>=0.5);
        else
            s0 = @(x,t) max(min(x,0.5) - max(x-t,0),0);
            s1 = @(x,t) max(x - max(x-t,0.5),0);
            exact = @(x,t) source_value*s0(x,t).*exp(-sigma_value*s1(x,t));
        end
    case 4 % void
        mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 1.0;
        phys.mu     = 1.0;
        phys.sigma  = @(x,t) 0.0;
        phys.source = @(x,t) 0.0;
        phys.speed  = 1;
        
        IC_option = 0;
        source_is_time_dependent = false;
        exact_solution_known = true;
        exact = @(x,t) t >= x;
    case 5 % MMS-1
        mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 0.0;
        phys.mu     = 1.0;
        phys.sigma  = @(x,t) 1.0;
        phys.speed  = 1;
        IC_option = 0;
        exact_solution_known = true;
        if temporal_scheme == 0 % steady-state
          source_is_time_dependent = false;
          phys.source = @(x,t) pi*cos(pi*x)+sin(pi*x);
          exact = @(x,t) sin(pi*x);
        else
          source_is_time_dependent = true;
          phys.source = @(x,t) sin(pi*x)+pi*t*cos(pi*x)+t*sin(pi*x);
          exact = @(x,t) t*sin(pi*x);
        end
    case 6 % MMS: TR: u = x*t, SS: u = x
        mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 0.0;
        phys.mu     = 1.0;
        phys.sigma  = @(x,t) 1.0;
        phys.speed  = 1;
        IC_option = 0;
        exact_solution_known = true;
        if temporal_scheme == 0 % steady-state
          source_is_time_dependent = false;
          phys.source = @(x,t) 1 + x;
          exact = @(x,t) x;
        else
          source_is_time_dependent = true;
          phys.source = @(x,t) x + t + x*t;
          exact = @(x,t) x*t;
        end
    otherwise
        error('Invalid problem ID chosen');
end

% set initial conditions function and compute initial solution
if ~(temporal_scheme == 0)
    switch IC_option
        case 0 % zero
            phys.IC = @(x) zeros(length(x),1);
        case 1 % exponential pulse
            phys.IC = @(x) exp(-200*(x-0.5).^2);
        case 2 % exponential and square pulse
            phys.IC = @(x) exp(-200*(x-0.3).^2)+(x>=0.6).*(x<=0.8);
        case 3 % square pulse
            phys.IC = @(x) (x>=0.4).*(x<0.7);
        otherwise
            error('Invalid IC option');
    end
end

% assert that periodic BC are not specified for a steady-state problem
assert(~(temporal_scheme == 0 & phys.periodic_BC),...
    'Cannot specify periodic BC for a steady-state problem.');

% assert that time-dependent source is not used for a steady-state problem
assert(~(temporal_scheme == 0 & source_is_time_dependent),...
    'Cannot have a time-dependent source in a steady-state problem.');

%% Setup

% determine if linear system needs to be modified for Dirichlet BC
modify_for_weak_DirichletBC   = ...
    ~phys.periodic_BC && ~impose_DirichletBC_strongly;
modify_for_strong_DirichletBC = ...
    ~phys.periodic_BC &&  impose_DirichletBC_strongly;

% group options
numerics_opts.modify_for_weak_DirichletBC   = modify_for_weak_DirichletBC;
numerics_opts.modify_for_strong_DirichletBC = modify_for_strong_DirichletBC;
numerics_opts.high_order_scheme             = high_order_scheme;
numerics_opts.theta                     = theta;
nonlin_solver_opts.tol                  = nonlin_tol;
nonlin_solver_opts.max_iter             = max_iter;
nonlin_solver_opts.relaxation_parameter = relaxation_parameter;

% number of dofs
if phys.periodic_BC
    dof_handler.n_dof = mesh.n_cell;
else
    dof_handler.n_dof = mesh.n_cell + 1;
end

% mesh
mesh.x = linspace(mesh.x_min,mesh.x_max,mesh.n_cell+1)'; % positions of dofs/cell vertices
mesh.dx = diff(mesh.x);                     % element sizes
if phys.periodic_BC
    mesh.x(end)=[];
end
x_center = 0.5*(mesh.x(1:end-1) + mesh.x(2:end));

% get quadrature points and weights and evaluate basis functions
[quadrature.zq,quadrature.wq]  = get_GL_quadrature(quadrature.nq);
[quadrature.v,quadrature.dvdz] = get_lagrange_basis(quadrature.zq,2); 
quadrature.Jac = 0.5*mesh.dx; % Jacobians of reference cell transformations

%% Material properties

% assert that a time-dependent source is not used; hasn't yet been
% implemented
assert(~source_is_time_dependent,...
    'Time-dependent source not yet implemented for DMP.');

% compute min and max sigma and source in the support of i for time 0
t = 0;
[sigma_min, sigma_max]  = compute_min_max_per_dof(...
    phys.sigma, t,dof_handler.n_dof,mesh,quadrature.zq);
[source_min,source_max] = compute_min_max_per_dof(...
    phys.source,t,dof_handler.n_dof,mesh,quadrature.zq);

%% Assembly

% create connectivity array
dof_handler.connectivity = [linspace(1,mesh.n_cell,mesh.n_cell)',...
    linspace(2,mesh.n_cell+1,mesh.n_cell)'];
% for periodic BC, last and first dofs are the same
if phys.periodic_BC
    dof_handler.connectivity(end,:) = [dof_handler.connectivity(end,1) 1];
end
% dofs per cell
dof_handler.dofs_per_cell = 2;
% number of cells
dof_handler.n_cell = mesh.n_cell;

% assemble mass matrices, inviscid and low-order steady-state matrices,
% low-order viscosity and corresponding artificial diffusion matrix
[MC,ML,A,AL,DL,viscL] = build_matrices(...
    quadrature,mesh,dof_handler,phys,...
    low_order_scheme,modify_for_weak_DirichletBC);

% assemble steady-state rhs at time 0
t = 0;
b = assemble_ss_rhs(t,quadrature,mesh,dof_handler,phys,...
    modify_for_weak_DirichletBC);

% check that speed is nonzero; otherwise, infinite dt will be computed
assert(phys.speed ~= 0,'Speed cannot be zero.');

% compute dt using CFL condition, even for implicit, just so that time
% steps will be equal between explicit and implicit, for comparison
max_speed_dx = 0.0;
for i = 1:dof_handler.n_dof
    max_speed_dx = max(max_speed_dx, abs(AL(i,i))/ML(i,i));
end
dtCFL = CFL / max_speed_dx;

% determine sytem matrices
if (temporal_scheme == 0) % steady-state
    system_matrix = AL;
else                      % transient
    system_matrix = ML;
end
% modify for Dirichlet BC
if modify_for_strong_DirichletBC
    system_matrix(1,:)=0; system_matrix(1,1)=1;
end

% create modified low-order ss matrix and ss rhs
AL_mod = AL;
b_mod = b;
if modify_for_strong_DirichletBC
    AL_mod(1,:)=0; AL_mod(1,1)=1;
    b_mod(1) = phys.inc;
end

%% Low-order Solution

if (compute_low_order)
    fprintf('\nComputing low-order solution...\n\n');
    
    if (temporal_scheme == 0) % steady-state
        % solve steady-state problem
        uL = AL_mod \ b_mod;
    else % transient
        % compute initial conditions
        u_old = phys.IC(mesh.x);
        if impose_BC_on_IC
          u_old(1) = phys.inc;
        end
        b_old = b;
        b_new = b_old;
        t = 0;
        time_step = 0;
        reached_end_of_transient = false;
        while (~reached_end_of_transient)
            % advance time step index
            time_step = time_step+1;
            
            % shorten time step if necessary
            if (t+dtCFL >= t_end)
                dt = t_end - t;
                reached_end_of_transient = true;
            else
                dt = dtCFL;
                reached_end_of_transient = false;
            end
            fprintf('Time step %i: t = %f->%f\n',time_step,t,t+dt);
            
            % perform step
            switch temporal_scheme
                case 0 % steady-state
                    error('This should not be reachable');
                case 1 % explicit Euler
                    % compute ss rhs
                    b_old = assemble_ss_rhs(t,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    
                    % compute solution
                    uL = compute_low_order_solution_theta(u_old,AL,...
                        ML,b_old,dt,0,phys.inc,modify_for_strong_DirichletBC);
                case 2 % SSP3
                    % stage 1
                    u_old_stage = u_old;
                    t_stage = t;
                    b_stage = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    uL_stage = compute_low_order_solution_theta(u_old_stage,AL,...
                        ML,b_stage,dt,0,phys.inc,...
                        modify_for_strong_DirichletBC);
                    % stage 2
                    u_old_stage = uL_stage;
                    t_stage = t + dt;
                    b_stage = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    uL_stage = compute_low_order_solution_theta(u_old_stage,AL,...
                        ML,b_stage,dt,0,phys.inc,...
                        modify_for_strong_DirichletBC);
                    % stage 3
                    u_old_stage = 0.75*u_old + 0.25*uL_stage;
                    t_stage = t + 0.5*dt;
                    b_stage = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    uL_stage = compute_low_order_solution_theta(u_old_stage,AL,...
                        ML,b_stage,dt,0,phys.inc,...
                        modify_for_strong_DirichletBC);
                    % final combination
                    uL = 1/3*u_old + 2/3*uL_stage;
                case 3 % theta
                    % compute ss rhs
                    b_new = assemble_ss_rhs(t+dt,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    b_theta = (1-theta)*b_old + theta*b_new;
                    
                    % compute solution
                    uL = compute_low_order_solution_theta(u_old,AL,...
                        ML,b_theta,dt,theta,phys.inc,...
                        modify_for_strong_DirichletBC);
                otherwise
                    error('Invalid temporal discretization scheme');
            end
            
            % test steady-state convergence
            ss_err = norm(uL-u_old,2);
            reached_steady_state = ss_err < ss_tol;
            if (reached_steady_state)
                fprintf('Transient terminated due to steady-state\n');
                break;
            end
            
            % plot
            if (plot_low_order_transient)
                plot(mesh.x,uL);
                legend('Low-order','Location',legend_location);
                axis([mesh.x_min mesh.x_max 0 1]);
                pause(pausetime);
            end
            
            % reset old quantities
            u_old = uL;
            b_old = b_new;
            t = t+dt;
        end % end transient loop
    end % end "is steady-state?" block
    
    % save uL for plotting because "uL" is used in FCT loop
    uL_final = uL;
end

%% High-order solution

if (compute_high_order)
    fprintf('\nComputing high-order solution...\n\n');
    
    if (temporal_scheme == 0) % steady-state
        [uH,DH,viscE] = compute_high_order_solution_ss(A,b,viscL,mesh,...
            phys,quadrature,ev,dof_handler,high_order_scheme,max_iter,...
            nonlin_tol,relaxation_parameter,modify_for_strong_DirichletBC);
    else % transient
        % compute initial conditions
        u_old = phys.IC(mesh.x);
        if impose_BC_on_IC
          u_old(1) = phys.inc;
        end
        u_older = u_old;
        b_old = b;
        b_new = b_old;
        t = 0;
        time_step = 0;
        dt_old = dtCFL;
        [DH,viscE] = compute_high_order_diffusion_matrix(u_old,...
            u_old,dt_old,viscL,mesh,phys,quadrature,ev,...
            dof_handler,high_order_scheme);
        AH = A + DH;
        
        % start transient loop
        reached_end_of_transient = false;
        while (~reached_end_of_transient)
            % advance time step index
            time_step = time_step+1;
            
            % shorten time step if necessary
            if (t+dtCFL >= t_end)
                dt = t_end - t;
                reached_end_of_transient = true;
            else
                dt = dtCFL;
                reached_end_of_transient = false;
            end
            fprintf('Time step %i: t = %f->%f\n',time_step,t,t+dt);
            
            switch temporal_scheme
                case 0 % steady-state
                    error('This should not be reachable');
                case 1 % explicit Euler
                    b = assemble_ss_rhs(t,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    [uH,DH] = high_order_step(u_older,u_old,dt_old,dt,...
                        A,b,MC,0,viscL,mesh,phys,quadrature,dof_handler,...
                        ev,high_order_scheme,...
                        modify_for_strong_DirichletBC);
                case 2 % SSP3
                    % stage 1
                    u_old_stage = u_old;
                    u_older_stage = u_older;
                    t_stage = t;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,...
                        high_order_scheme,modify_for_strong_DirichletBC);
                    % stage 2
                    u_old_stage = uH_stage;
                    u_older_stage = u_older;
                    t_stage = t + dt;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,...
                        high_order_scheme,modify_for_strong_DirichletBC);
                    % stage 3
                    u_old_stage = 0.75*u_old + 0.25*uH_stage;
                    u_older_stage = u_older;
                    t_stage = t + 0.5*dt;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,...
                        high_order_scheme,modify_for_strong_DirichletBC);
                    % final combination
                    uH = 1/3*u_old + 2/3*uH_stage;
                case 3 % theta
                    % compute new ss rhs
                    b_new = assemble_ss_rhs(t+dt,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);

                    % compute solution
                    [uH,AH,DH] = compute_high_order_solution_theta(...
                        u_old,dt,MC,A,AH,b_old,b_new,viscL,quadrature,...
                        mesh,dof_handler,phys,ev,numerics_opts,...
                        nonlin_solver_opts);
                otherwise
                    error('Invalid temporal discretization scheme');
            end

            % check for NaNs in solution
            check_NaN(uH);
            
            % test steady-state convergence
            ss_err = norm(uH-u_old,2);
            reached_steady_state = ss_err < ss_tol;
            if (reached_steady_state)
                fprintf('Transient terminated due to steady-state\n');
                break;
            end
            
            % plot
            if (plot_high_order_transient)
                plot(mesh.x,uH);
                axis([-inf inf -inf 0.5]);
                legend('High-order');
                pause(pausetime);
            end
            
            % reset u_old and u_older and advance time
            u_older = u_old;
            u_old   = uH;
            b_old   = b_new;
            dt_old = dt;
            t = t+dt;
        end % end transient loop
    end % end "is steady-state?" block
    
    % save uH for plotting because "uH" is used in FCT loop
    uH_final = uH;
end

%% FCT solution

if (compute_FCT)
    fprintf('\nComputing FCT solution...\n\n');
    
    if (temporal_scheme == 0) % steady-state
        % compute high-order ss solution
        [uH,DH,viscE] = compute_high_order_solution_ss(A,b,viscL,mesh,...
            phys,quadrature,ev,dof_handler,high_order_scheme,max_iter,...
            nonlin_tol,relaxation_parameter,modify_for_strong_DirichletBC);
        
        % compute flux correction matrix
        F = flux_correction_matrix_ss(uH,DL-DH);
        fluxvector = convert_matrix_to_edge_vector(F);
        
        % compute low-order solution
        uL = AL_mod \ b_mod;

        % prelimit flux corrections if user specified
        if (prelimit)
            F = prelimit_fluxes(F,uL);
        end
        
        % initialize solution iterate
        switch FCT_initialization
            case 1 % zero
                uFCT = zeros(dof_handler.n_dof,1);
            case 2 % low-order
                uFCT = uL;
            case 3 % high-order
                uFCT = uH;
            otherwise
                error('Invalid FCT initialization option');
        end
        
        % iteration loop
        for iter = 1:max_iter
            % compute limited flux correction sum
            [flim,Wminus,Wplus] = compute_limited_flux_sums_ss(uFCT,F,...
                AL_mod,b_mod,...
                sigma_min,sigma_max,source_min,source_max,phys,...
                dof_handler.n_dof,DMP_option,limiting_option);
            
            % compute correction rhs and modify for Dirichlet BC
            system_rhs = b + flim;
            if (impose_DirichletBC_strongly)
                system_rhs(1) = phys.inc;
            end
            
            % compute residual
            ss_res = system_rhs - AL_mod*uFCT;
            
            % test convergence of previous iteration
            converged = test_convergence(ss_res,nonlin_tol,iter);
            if (converged)
                fprintf('\tConverged at iteration %i\n',iter);
                break;
            end

            % solve modified system
            du = AL_mod \ ss_res;
            uFCT = uFCT + du;

            % plot
            if (plot_FCT_transient)
                figure(1);
                clf;
                hold on;
                plot(mesh.x,uH,'rx');
                plot(mesh.x,uL,'b+');
                plot(mesh.x,uFCT,'g');
                plot(mesh.x,Wminus,'k.');
                plot(mesh.x,Wplus,'ko');
                legend('High','Low','FCT(l+1)','W-(l)','W+(l)',...
                    'Location',legend_location);
                w = waitforbuttonpress;
            end
        end
        
        if (~converged)
            error('FCT Solution did not converge');
        end
    else
        % compute initial conditions
        u_old = phys.IC(mesh.x);
        if impose_BC_on_IC
          u_old(1) = phys.inc;
        end
        u_older = u_old;
        b_old = b;
        b_new = b_old;
        t = 0;
        time_step = 0;
        dt_old = dtCFL;
        [DH,viscE] = compute_high_order_diffusion_matrix(u_old,...
            u_old,dt_old,viscL,mesh,phys,quadrature,ev,...
            dof_handler,high_order_scheme);
        AH = A + DH;
        
        reached_end_of_transient = false;
        while (~reached_end_of_transient)
            % advance time step index
            time_step = time_step+1;
            
            % shorten time step if necessary
            if (t+dtCFL >= t_end)
                dt = t_end - t;
                reached_end_of_transient = true;
            else
                dt = dtCFL;
                reached_end_of_transient = false;
            end
            fprintf('Time step %i: t = %f->%f\n',time_step,t,t+dt);
            
            switch temporal_scheme
                case 0 % steady-state
                    error('This should not be reachable');
                case 1 % explicit Euler
                    % perform high-order step
                    b = assemble_ss_rhs(t,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    [uH,DH] = high_order_step(u_older,u_old,dt_old,dt,...
                        A,b,MC,0,viscL,mesh,phys,quadrature,dof_handler,...
                        ev,high_order_scheme,...
                        modify_for_strong_DirichletBC);
                    
                    % perform FCT step
                    uFCT = FCT_step_explicit(u_old,uH,dt,...
                        ML,MC,AL,DH,DL,b,phys.inc,phys.speed,sigma_min,sigma_max,...
                        source_min,source_max,DMP_option,limiting_option,...
                        phys.periodic_BC,modify_for_strong_DirichletBC,prelimit);
                case 2 % SSP3
                    % stage 1
                    u_old_stage = u_old;
                    u_older_stage = u_older;
                    t_stage = t;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,...
                        high_order_scheme,modify_for_strong_DirichletBC);
                    uFCT_stage = FCT_step_explicit(u_old_stage,uH_stage,dt,...
                        ML,MC,AL,DH,DL,b,phys.inc,phys.speed,sigma_min,sigma_max,...
                        source_min,source_max,DMP_option,limiting_option,...
                        phys.periodic_BC,modify_for_strong_DirichletBC,prelimit);
                    
                   % stage 2
                    u_old_stage = uFCT_stage;
                    u_older_stage = u_older;
                    t_stage = t + dt;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,...
                        high_order_scheme,modify_for_strong_DirichletBC);
                    uFCT_stage = FCT_step_explicit(u_old_stage,uH_stage,dt,...
                        ML,MC,AL,DH,DL,b,phys.inc,phys.speed,sigma_min,sigma_max,...
                        source_min,source_max,DMP_option,limiting_option,...
                        phys.periodic_BC,modify_for_strong_DirichletBC,prelimit);
                    
                    % stage 3
                    u_old_stage = 0.75*u_old + 0.25*uFCT_stage;
                    u_older_stage = u_older;
                    t_stage = t + 0.5*dt;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,...
                        high_order_scheme,modify_for_strong_DirichletBC);
                    uFCT_stage = FCT_step_explicit(u_old_stage,uH_stage,dt,...
                        ML,MC,AL,DH,DL,b,phys.inc,phys.speed,sigma_min,sigma_max,...
                        source_min,source_max,DMP_option,limiting_option,...
                        phys.periodic_BC,modify_for_strong_DirichletBC,prelimit);
                    
                    % final combination
                    uFCT = 1/3*u_old + 2/3*uFCT_stage;
                case 3 % theta
                    % compute old steady-state residual
                    ss_res = b_old - AL*u_old;
                    
                    % compute system matrix
                    system_matrix = ML + theta*dt*AL;
                    if (modify_for_strong_DirichletBC)
                        system_matrix(1,:)=0; system_matrix(1,1)=1;
                    end
                    
                    % compute new ss rhs
                    b_new = assemble_ss_rhs(t+dt,quadrature,mesh,...
                        dof_handler,phys,modify_for_weak_DirichletBC);
                    
                    % compute high-order solution
                    [uH,AH,DH] = compute_high_order_solution_theta(...
                        u_old,dt,MC,A,AH,b_old,b_new,viscL,quadrature,...
                        mesh,dof_handler,phys,ev,numerics_opts,...
                        nonlin_solver_opts);

                    % compute flux correction matrix
                    F = flux_correction_matrix(u_old,uH,dt,DH,DL,MC,theta);
                    
                    % compute theta ss rhs
                    b_theta = (1-theta)*b_old + theta*b_new;
                    
                    % compute solution
                    uL = compute_low_order_solution_theta(u_old,AL,...
                        ML,b_theta,dt,theta,phys.inc,...
                        modify_for_strong_DirichletBC);
                    
                    % initialize solution iterate
                    %uFCT = uL;
                    uFCT = zeros(dof_handler.n_dof,1);
                    
                    % iteration loop
                    for iter = 1:max_iter
                        % compute limited flux correction sum
                        flim = compute_limited_flux_sums(u_old,uFCT,dt,...
                            ML,AL,b,F,sigma_min,sigma_max,source_min,...
                            source_max,theta,dof_handler.n_dof,...
                            phys.speed,phys.inc,phys.periodic_BC,...
                            limiting_option,DMP_option);
                        
                        % compute system rhs
                        system_rhs = ML*u_old + (1-theta)*dt*ss_res ...
                            + theta*dt*b_new + dt*flim;
                        if (modify_for_strong_DirichletBC)
                            system_rhs(1) = phys.inc;
                        end
                        
                        % compute residual
                        res = system_rhs - system_matrix*uFCT;
                        
                        % test convergence of previous iteration
                        converged = test_convergence(res,nonlin_tol,iter);
                        if (converged)
                            fprintf('\tConverged at iteration %i\n',iter);
                            break;
                        end
                        
                        % compute change in solution iterate
                        du = system_matrix \ system_rhs - uFCT;
    
                        % solve modified system
                        uFCT = uFCT + relaxation_parameter * du;
                    end
                    
                    % report if the solution did not converge
                    if (~converged)
                        error('FCT Solution did not converge');
                    end
                otherwise
                    error('Invalid temporal discretization scheme');
            end
            
            % test steady-state convergence
            ss_err = norm(uFCT-u_old,2);
            reached_steady_state = ss_err < ss_tol;
            if (reached_steady_state)
                fprintf('Transient terminated due to steady-state\n');
                break;
            end
            
            % plot
            if (plot_FCT_transient)
                plot(mesh.x,uFCT,'g');
                legend('FCT','Location',legend_location);
                pause(pausetime);
            end
            
            % reset u_old and u_older and advance time
            u_older = u_old;
            u_old   = uFCT;
            b_old   = b_new;
            dt_old = dt;
            t = t+dt;
        end
    end
end

%% Plot Solution

figure; clf; hold on;

% plot exact solution
if (exact_solution_known)
    xx = linspace(mesh.x_min,mesh.x_max,1000)';
    u_exact = exact(xx,t_end);
    plot(xx,u_exact,'k-');
    legend_entries = char('Exact');
end

% plot low-order solution
if (compute_low_order)
    plot(mesh.x,uL_final,'r-s');
    if (exact_solution_known)
        legend_entries = char(legend_entries,'Low-order');
    else
        legend_entries = char('Low-order');
    end
end

% plot high-order solution
if (compute_high_order)
    plot(mesh.x,uH_final,'b-+');
    legend_entries = char(legend_entries,'High-order');
end

% plot FCT solution
if (compute_FCT)
    plot(mesh.x,uFCT,'g-x');
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
legend(legend_entries,'Location',legend_location);

%% Plot viscosity

% plot viscosity if requested
if (plot_viscosity)
    figure; clf;
    
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
    legend(legend_entries,'Location',legend_location);
end

%% Output

% save exact solution
if (save_exact_solution)
    exact_file = 'output/uexact.txt';
    dlmwrite(exact_file,[xx,u_exact],' ');
end

% determine string to be appended to results for time discretization
switch temporal_scheme
    case 0
        time_string = 'ss';
    case 1
        time_string = 'FE';
    case 2
        time_string = 'SSPRK33';
    case 3
        time_string = sprintf('theta%0.1f',theta);
        % rename common theta schemes
        if (strcmp(time_string,'theta0.0'))
            time_string = 'FE';
        elseif (strcmp(time_string,'theta1.0'))
            time_string = 'BE';
        elseif (strcmp(time_string,'theta0.5'))
            time_string = 'CN';
        end
    otherwise
        error('Invalid time discretization scheme');
end

% save low-order solution
if (save_low_order_solution)
    if (compute_low_order)
        low_order_file = ['output/uL_',time_string,'.txt'];
        dlmwrite(low_order_file,[mesh.x,uL_final],' ');
    end
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
    if (compute_high_order)
        high_order_file = ['output/uH_',high_order_string,'_',time_string,'.txt'];
        dlmwrite(high_order_file,[mesh.x,uH_final],' ');
    end
end

% save FCT solution
if (save_FCT_solution)
    if (compute_FCT)
        FCT_file = ['output/uFCT_',high_order_string,'_',time_string,'.txt'];
        dlmwrite(FCT_file,[mesh.x,uFCT],' ');
    end
end
