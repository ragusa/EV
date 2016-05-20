% optional parameter: number of cells
function [return_value1,return_value2] = main(n_cell)

clc; close all; 
if (nargin == 0)
  clear;
end

%% User Options
%--------------------------------------------------------------------------
% finite element options
%--------------------------------------------------------------------------
if (nargin > 0)
  mesh.n_cell = n_cell; % number of elements
else
  mesh.n_cell = 32;     % number of elements
end
quadrature.nq = 3;      % number of quadrature points per cell
opts.impose_DirichletBC_strongly = false; % impose Dirichlet BC strongly?
opts.use_penalty_bc = true; % add penalty BC?
opts.penalty_bc_coef = 1000.0; % penalty BC coef: alpha*A(i,i) = alpha*inc
%--------------------------------------------------------------------------
% spatial method options
%--------------------------------------------------------------------------
compute_low_order  = true; % compute and plot low-order solution?
compute_high_order = true; % compute and plot high-order solution?
compute_FCT        = true; % compute and plot FCT solution?

% low_order_scheme: 1 = algebraic low-order scheme
%                   2 = graph-theoretic low-order scheme
opts.low_order_scheme  = 2;

% high_order_scheme: 1 = Galerkin
%                    2 = Entropy viscosity
%                    3 = Alternate Entropy viscosity 1
%                    4 = Alternate Entropy viscosity 2 (should be same as 1)
opts.high_order_scheme = 2;

%--------------------------------------------------------------------------
% entropy viscosity options
%--------------------------------------------------------------------------
ev.cE = 0.1; % coefficient for entropy residual in entropy viscosity
ev.cJ = ev.cE*1; % coefficient for jumps in entropy viscosity
ev.entropy       = @(u) 0.5*u.^2; % entropy function
ev.entropy_deriv = @(u) u;        % derivative of entropy function
ev.use_local_ev_norm = false; % option to use local entropy normalization
ev.smooth_entropy_viscosity = false; % option to smooth entropy viscosity
ev.smoothing_weight = 1.0; % weight for center value in smoothing
%--------------------------------------------------------------------------
% time options
%--------------------------------------------------------------------------
% temporal scheme: 0 = steady-state
%                  1 = SSPRK(1,1) (Explicit Euler)
%                  2 = SSPRK(3,3) (Shu-Osher)
%                  3 = theta method
opts.temporal_scheme = 0; % temporal discretization scheme

% theta parameter to use if using a theta method: 0.0 = FE
%                                                 0.5 = CN
%                                                 1.0 = BE
opts.theta = 1.0;     
       
opts.use_constant_dt = false; % option to use constant dt instead of CFL
opts.constant_dt = 0.001;    % time step size to use if using constant size
opts.CFL = 10.0;       % CFL number
opts.t_end = 1.0;     % max time to run
opts.ss_tol = 1.0e-6;  % steady-state tolerance
%--------------------------------------------------------------------------
% FCT options
%--------------------------------------------------------------------------
% DMP option: 1 = low-order DMP
%             2 = widen low-order DMP to analytic
%             3 = analytic
%             4 = analytic upwind
%             5 = analytic alternate
fct_opts.DMP_option = 4;

% limiter option: 0 = All 0 (no correction; low-order)
%                 1 = All 1 (full correction; high-order)
%                 2 = Zalesak limiter
%                 3 = Josh limiter
fct_opts.limiting_option = 2;

fct_opts.use_multipass_limiting = false; % option to use multi-pass limiting
fct_opts.multipass_tol = 0.01; % percentage tolerance for multi-pass

% option to enforce Q+ >= 0, Q- <= 0
fct_opts.enforce_antidiffusion_bounds_signs = true;

% FCT initialization option: 1 = zeros
%                            2 = low-order solution
%                            3 = high-order solution
fct_opts.FCT_initialization = 2;

% option to skip limitation of bounds if solution bounds are satisfied already
fct_opts.skip_limiter_if_bounds_satisfied = false;

% option to prelimit correction fluxes
fct_opts.prelimit = false;

% limiting coefficient bounds for Dirichlet nodes
fct_opts.dirichlet_limiting_coefficient = 0.0; 
%--------------------------------------------------------------------------
% physics options
%--------------------------------------------------------------------------
% problemID: 0: custom - use parameters below
%            1: pure absorber without source
%            2: void without source -> absorber without source
%            3: void with    source -> absorber without source
%            4: void
%            5: MMS: TR: u = t*sin(pi*x)  SS: u = sin(pi*x)
%            6: MMS: TR: u = x*t          SS: u = x
%            7: source_in_absorber
%            8: interface
problemID = 2;

% IC_option: 0: zero
%            1: exponential pulse
%            2: exponential and square pulse
IC_option = 0;

mesh.x_min = 0.0;         % left end of domain
mesh.x_max = 1.0;         % right end of domain

phys.periodic_BC = false; % option for periodic BC; otherwise Dirichlet
phys.mu     = 1;          % cos(angle)
sigma_value  = 100.0;
source_value = 10.0;
phys.sigma  = @(x,t) sigma_value; % cross section function
phys.source = @(x,t) (x<0.5)*source_value; % source function
phys.inc    = 0;          % incoming flux
phys.speed  = 1;          % advection speed
phys.source_is_time_dependent = false; % is source time-dependent?
phys.impose_BC_on_IC = true; % option to impose Dirichlet BC on IC
exact_solution_known = true;
exact = @(x,t) (x<0.5).*(source_value/sigma_value*(1.0-exp(-sigma_value*x)))...
  + (x>=0.5).*(source_value/sigma_value*(1.0-exp(-sigma_value*0.5))...
  * exp(-sigma_value*(x-0.5)));
%--------------------------------------------------------------------------
% nonlinear solver options
%--------------------------------------------------------------------------
nonlin_opts.max_iter = 10000;    % maximum number of nonlinear solver iterations
nonlin_opts.nonlin_tol = 1e-10; % nonlinear solver tolerance for discrete L2 norm
nonlin_opts.relax = 1.0; % relaxation parameter for iteration
%--------------------------------------------------------------------------
% plot options
%--------------------------------------------------------------------------
out_opts.plot_low_order_transient  = false; % plot low-order transient?
out_opts.plot_high_order_transient = false; % plot high-order transient?
out_opts.plot_FCT_transient        = false; % plot FCT transient?

out_opts.plot_EV_iteration         = false; % plot EV iteration?
out_opts.plot_FCT_iteration        = false; % plot FCT iteration?

out_opts.plot_viscosity            = false; % plot viscosities?

out_opts.pause_type                = 'wait'; % pause type: 'wait' or 'time'
out_opts.pausetime                 = 0.5; % time to pause for transient plots
out_opts.legend_location           = 'NorthEast'; % location of plot legend
%--------------------------------------------------------------------------
% output options
%--------------------------------------------------------------------------
% option to output l2 norm of entropy residual
return_value_option = 0; % 0: nothing - just return zero
                         % 1: L^2 norm of entropy residual
                         % 2: L^2 norm of entropy jumps

save_exact_solution      = true; % option to save exact solution 
save_low_order_solution  = true; % option to save low-order solution
save_high_order_solution = true; % option to save high-order solution
save_FCT_solution        = true; % option to save FCT solution
save_FCT_bounds          = true; % option to save FCT bounds
save_antidiffusion_matrix = false; % option to save antidiffusion matrix
%-------------------------------------------------------------------------

%% Define Problem

% set flag for steady state
opts.is_steady_state = opts.temporal_scheme == 0;

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
        sigma_value = 10.0;
        phys.sigma  = @(x,t) sigma_value;
        phys.source = @(x,t) 0.0;
        phys.speed  = 1;
        
        IC_option = 0;
        phys.source_is_time_dependent = false;
        exact_solution_known = true;
        if opts.is_steady_state
          exact = @(x,t) exp(-sigma_value*x);
        else
          exact = @(x,t) (t>=x).*exp(-sigma_value*x);
        end
    case 2 % void without source -> absorber without source
        mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 1.0;
        phys.mu     = 1.0;
        sigma_value = 100.0;
        phys.sigma  = @(x,t) sigma_value*(x >= 0.5);
        phys.source = @(x,t) 0.0;
        phys.speed  = 1;
        
        IC_option = 0;
        phys.source_is_time_dependent = false;
        exact_solution_known = true;
        if opts.is_steady_state
          exact = @(x,t) (x<0.5) + (x>=0.5).*(exp(-sigma_value*(x-0.5)));
        else
          exact = @(x,t) (t>=x).*((x<0.5) + (x>=0.5).*(exp(-sigma_value*(x-0.5))));
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
        phys.source_is_time_dependent = false;
        exact_solution_known = true;
        if opts.is_steady_state
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
        phys.source_is_time_dependent = false;
        exact_solution_known = true;
        exact = @(x,t) t >= x;
    case 5 % MMS: TR: u = t*sin(pi*x), SS: u = sin(pi*x)
        mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 0.0;
        phys.mu     = 1.0;
        phys.sigma  = @(x,t) 1.0;
        phys.speed  = 1;
        IC_option = 0;
        exact_solution_known = true;
        if opts.is_steady_state % steady-state
          phys.source_is_time_dependent = false;
          phys.source = @(x,t) pi*cos(pi*x)+sin(pi*x);
          exact = @(x,t) sin(pi*x);
        else
          phys.source_is_time_dependent = true;
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
        if opts.is_steady_state % steady-state
          phys.source_is_time_dependent = false;
          phys.source = @(x,t) 1 + x;
          exact = @(x,t) x;
        else
          phys.source_is_time_dependent = true;
          phys.source = @(x,t) x + t + x*t;
          exact = @(x,t) x*t;
        end
    case 7 % source_in_absorber
         mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 0;
        phys.mu     = 1;
        source_value = 10.0;
        sigma_value  = 100.0;
        phys.source = @(x,t) (x<0.5)*source_value;
        phys.sigma  = @(x,t) sigma_value;
        phys.speed  = 1;
        IC_option = 0;
        phys.source_is_time_dependent = false;
        exact_solution_known = true;
        if opts.is_steady_state % steady-state
            exact = @(x,t) (x<0.5).*(source_value/sigma_value*(1.0-exp(-sigma_value*x)))...
                + (x>=0.5).*(source_value/sigma_value*(1.0-exp(-sigma_value*0.5))...
                * exp(-sigma_value*(x-0.5)));
        else
            error('Exact solution not written for transient case');
        end
    case 8 % interface
        mesh.x_min = 0.0;
        mesh.x_max = 1.0;
        phys.periodic_BC = false;
        phys.inc    = 0;
        phys.mu     = 1;
        x1 = 0.5;
        source_value1 = 10.0;
        source_value2 = 20.0;
        sigma_value1  = 10.0;
        sigma_value2  = 40.0;
        phys.source = @(x,t) (x<x1)*source_value1 + (x>=x1)*source_value2;
        phys.sigma  = @(x,t) (x<x1)*sigma_value1  + (x>=x1)*sigma_value2;
        phys.speed  = 1;
        IC_option = 0;
        phys.source_is_time_dependent = false;
        exact_solution_known = true;
        if opts.is_steady_state % steady-state
            exact = @(x,t) (x<x1).*(source_value1/sigma_value1*(1.0-exp(-sigma_value1*x)))...
                + (x>=x1).*((source_value1/sigma_value1*(1.0-exp(-sigma_value1*x1))...
                * exp(-sigma_value2*(x-x1)))...
                + source_value2/sigma_value2*(1.0-exp(-sigma_value2*(x-x1))));
        else
            error('Exact solution not written for transient case');
        end
    otherwise
        error('Invalid problem ID chosen');
end

% set initial conditions function and compute initial solution
if ~(opts.is_steady_state)
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
assert(~(opts.is_steady_state & phys.periodic_BC),...
    'Cannot specify periodic BC for a steady-state problem.');

% assert that time-dependent source is not used for a steady-state problem
assert(~(opts.is_steady_state & phys.source_is_time_dependent),...
    'Cannot have a time-dependent source in a steady-state problem.');

% compute limiting coefficients
switch fct_opts.limiting_option
    case 0 % Full limiter
        fct_opts.limiter = @limiter_zeroes;
    case 1 % No limiter
        fct_opts.limiter = @limiter_ones;
    case 2 % Zalesak limiter
        fct_opts.limiter = @limiter_zalesak;
    case 3 % Josh limiter
        fct_opts.limiter = @limiter_josh;
    otherwise
        error('Invalid limiting option');
end

%% Setup

% determine if linear system needs to be modified for Dirichlet BC
opts.modify_for_weak_DirichletBC   = ...
    ~phys.periodic_BC && ~opts.impose_DirichletBC_strongly;
opts.modify_for_strong_DirichletBC = ...
    ~phys.periodic_BC &&  opts.impose_DirichletBC_strongly;

% number of dofs
if phys.periodic_BC
    dof_handler.n_dof = mesh.n_cell;
else
    dof_handler.n_dof = mesh.n_cell + 1;
end

% mesh
% positions of dofs/cell vertices
mesh.x = linspace(mesh.x_min,mesh.x_max,mesh.n_cell+1)'; 
if phys.periodic_BC
    mesh.x(end)=[];
end
% cell center positions
mesh.x_center = 0.5*(mesh.x(1:end-1) + mesh.x(2:end));
% element sizes
mesh.dx = diff(mesh.x);
% minimum distance
mesh.dx_min = min(mesh.dx);

% get quadrature points and weights and evaluate basis functions
[quadrature.zq,quadrature.wq]  = get_GL_quadrature(quadrature.nq);
[quadrature.v,quadrature.dvdz] = get_lagrange_basis(quadrature.zq,2); 
quadrature.Jac = 0.5*mesh.dx; % Jacobians of reference cell transformations

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
    quadrature,mesh,dof_handler,phys,opts);

% assemble steady-state rhs at time 0
t = 0;
b = assemble_ss_rhs(t,quadrature,mesh,dof_handler,phys,...
    opts.modify_for_weak_DirichletBC);

% modify steady-state RHS
if (opts.use_penalty_bc)
    % currently this is only implemented for steady-state; if not steady-
    % state, then one will need to perform this modification on the total
    % system RHS, not just the steady-state RHS
    if (opts.temporal_scheme ~= 0)
        error('Penalty BC implemented only for steady-state');
    end
    
    % modify RHS
    b(1) = b(1) + opts.penalty_bc_coef * phys.inc;
end

% check that speed is nonzero; otherwise, infinite dt will be computed
assert(phys.speed ~= 0,'Speed cannot be zero.');

% choose nominal time step size (last time step may be shortened)
if (opts.use_constant_dt) % use constant dt
    dt_nominal = opts.constant_dt;
else % compute dt from CFL
    % compute dt using CFL condition, even for implicit, just so that time
    % steps will be equal between explicit and implicit, for comparison
    max_speed_dx = 0.0;
    for i = 1:dof_handler.n_dof
        max_speed_dx = max(max_speed_dx, abs(AL(i,i))/ML(i,i));
    end
    dt_nominal = opts.CFL / max_speed_dx;
end

% determine sytem matrices
if (opts.is_steady_state) % steady-state
    system_matrix = AL;
else                      % transient
    system_matrix = ML;
end
% modify for Dirichlet BC
if opts.modify_for_strong_DirichletBC
    system_matrix(1,:)=0; system_matrix(1,1)=1;
end

% create modified low-order ss matrix and ss rhs
AL_mod = AL;
b_mod = b;
if opts.modify_for_strong_DirichletBC
    AL_mod(1,:)=0; AL_mod(1,1)=1;
    b_mod(1) = phys.inc;
end

%% Material properties

% assert that a time-dependent source is not used; hasn't yet been
% implemented
assert(~phys.source_is_time_dependent,...
    'Time-dependent source not yet implemented for DMP.');

% compute max distance traveled in a time step. if steady-state, take
% some value between 0 and dx
distance = 0.0;
if (opts.temporal_scheme == 0)
    distance = mesh.dx(1)*0.5;
else
    distance = phys.speed*dt_nominal;
end

% create source/sigma function
source_over_sigma = @(x,t) phys.source(x,t) ./ phys.sigma(x,t);
    
% compute min and max sigma and source in the support of i for time 0
t = 0;
if (fct_opts.DMP_option == 4) % upwind maximum principle
    [sigma_min, sigma_max]  = compute_min_max_per_dof_upwind(...
        phys.sigma, t,dof_handler.n_dof,mesh,quadrature.zq,distance);
    [source_min,source_max] = compute_min_max_per_dof_upwind(...
        phys.source,t,dof_handler.n_dof,mesh,quadrature.zq,distance);
    [source_over_sigma_min, source_over_sigma_max] = ...
        compute_min_max_per_dof_upwind(...
        source_over_sigma,t,dof_handler.n_dof,mesh,quadrature.zq,distance);
else
    [sigma_min, sigma_max]  = compute_min_max_per_dof(...
        phys.sigma, t,dof_handler.n_dof,mesh,quadrature.zq,distance);
    [source_min,source_max] = compute_min_max_per_dof(...
        phys.source,t,dof_handler.n_dof,mesh,quadrature.zq,distance);
    [source_over_sigma_min, source_over_sigma_max] = ...
        compute_min_max_per_dof(...
        source_over_sigma,t,dof_handler.n_dof,mesh,quadrature.zq,distance);
end

%% Low-order Solution

if (compute_low_order)
    fprintf('\nComputing low-order solution...\n\n');
    
    if (opts.is_steady_state) % steady-state
        % solve steady-state problem
        uL = AL_mod \ b_mod;
    else % transient
        % compute initial conditions
        u_old = phys.IC(mesh.x);
        if phys.impose_BC_on_IC
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
            if (t+dt_nominal >= opts.t_end)
                dt = opts.t_end - t;
                reached_end_of_transient = true;
            else
                dt = dt_nominal;
                reached_end_of_transient = false;
            end
            fprintf('Time step %i: t = %f->%f\n',time_step,t,t+dt);
            
            % perform step
            switch opts.temporal_scheme
                case 0 % steady-state
                    error('This should not be reachable');
                case 1 % explicit Euler
                    % compute ss rhs
                    b_old = assemble_ss_rhs(t,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    
                    % compute solution
                    uL = compute_low_order_solution_theta(u_old,AL,...
                        ML,b_old,dt,0,phys.inc,opts.modify_for_strong_DirichletBC);
                case 2 % SSP3
                    % stage 1
                    u_old_stage = u_old;
                    t_stage = t;
                    b_stage = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    uL_stage = compute_low_order_solution_theta(u_old_stage,AL,...
                        ML,b_stage,dt,0,phys.inc,...
                        opts.modify_for_strong_DirichletBC);
                    % stage 2
                    u_old_stage = uL_stage;
                    t_stage = t + dt;
                    b_stage = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    uL_stage = compute_low_order_solution_theta(u_old_stage,AL,...
                        ML,b_stage,dt,0,phys.inc,...
                        opts.modify_for_strong_DirichletBC);
                    % stage 3
                    u_old_stage = 0.75*u_old + 0.25*uL_stage;
                    t_stage = t + 0.5*dt;
                    b_stage = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    uL_stage = compute_low_order_solution_theta(u_old_stage,AL,...
                        ML,b_stage,dt,0,phys.inc,...
                        opts.modify_for_strong_DirichletBC);
                    % final combination
                    uL = 1/3*u_old + 2/3*uL_stage;
                case 3 % theta
                    % compute ss rhs
                    b_new = assemble_ss_rhs(t+dt,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    b_theta = (1-opts.theta)*b_old + opts.theta*b_new;
                    
                    % compute solution
                    uL = compute_low_order_solution_theta(u_old,AL,...
                        ML,b_theta,dt,opts.theta,phys.inc,...
                        opts.modify_for_strong_DirichletBC);
                otherwise
                    error('Invalid temporal discretization scheme');
            end
            
            % test steady-state convergence
            ss_err = norm(uL-u_old,2);
            reached_steady_state = ss_err < opts.ss_tol;
            if (reached_steady_state)
                fprintf('Transient terminated due to steady-state\n');
                break;
            end
            
            % plot
            if (out_opts.plot_low_order_transient)
                plot(mesh.x,uL);
                legend('Low-order','Location',out_opts.legend_location);
                axis([mesh.x_min mesh.x_max 0 1]);
                pause(out_opts.pausetime);
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
    
    if (opts.is_steady_state) % steady-state
        [uH,DH,viscE] = compute_high_order_solution_ss(A,b,viscL,mesh,...
            phys,quadrature,ev,dof_handler,opts,nonlin_opts,out_opts);
    else % transient
        % compute initial conditions
        u_old = phys.IC(mesh.x);
        if phys.impose_BC_on_IC
          u_old(1) = phys.inc;
        end
        u_older = u_old;
        b_old = b;
        b_new = b_old;
        t = 0;
        time_step = 0;
        dt_old = dt_nominal;
        [DH,viscE] = compute_high_order_diffusion_matrix(u_old,...
            u_old,dt_old,viscL,mesh,phys,quadrature,ev,...
            dof_handler,opts.high_order_scheme);
        AH = A + DH;
        
        % start transient loop
        reached_end_of_transient = false;
        while (~reached_end_of_transient)
            % advance time step index
            time_step = time_step+1;
            
            % shorten time step if necessary
            if (t+dt_nominal >= opts.t_end)
                dt = opts.t_end - t;
                reached_end_of_transient = true;
            else
                dt = dt_nominal;
                reached_end_of_transient = false;
            end
            fprintf('Time step %i: t = %f->%f\n',time_step,t,t+dt);
            
            switch opts.temporal_scheme
                case 0 % steady-state
                    error('This should not be reachable');
                case 1 % explicit Euler
                    b = assemble_ss_rhs(t,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    [uH,DH] = high_order_step(u_older,u_old,dt_old,dt,...
                        A,b,MC,0,viscL,mesh,phys,quadrature,dof_handler,...
                        ev,opts);
                case 2 % SSP3
                    % stage 1
                    u_old_stage = u_old;
                    u_older_stage = u_older;
                    t_stage = t;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,opts);
                    % stage 2
                    u_old_stage = uH_stage;
                    u_older_stage = u_older;
                    t_stage = t + dt;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,opts);
                    % stage 3
                    u_old_stage = 0.75*u_old + 0.25*uH_stage;
                    u_older_stage = u_older;
                    t_stage = t + 0.5*dt;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,opts);
                    % final combination
                    uH = 1/3*u_old + 2/3*uH_stage;
                case 3 % theta
                    % compute new ss rhs
                    b_new = assemble_ss_rhs(t+dt,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);

                    % compute solution
                    [uH,AH,DH] = compute_high_order_solution_theta(...
                        u_old,dt,MC,A,AH,b_old,b_new,viscL,quadrature,...
                        mesh,dof_handler,phys,ev,opts,...
                        nonlin_opts);
                otherwise
                    error('Invalid temporal discretization scheme');
            end

            % check for NaNs in solution
            check_NaN(uH);
            
            % test steady-state convergence
            ss_err = norm(uH-u_old,2);
            reached_steady_state = ss_err < opts.ss_tol;
            if (reached_steady_state)
                fprintf('Transient terminated due to steady-state\n');
                break;
            end
            
            % plot
            if (out_opts.plot_high_order_transient)
                plot(mesh.x,uH);
                axis([-inf inf -inf 0.5]);
                legend('High-order');
                pause(out_opts.pausetime);
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
    
    if (opts.is_steady_state) % steady-state
        % compute high-order ss solution
        [uH,DH,viscE] = compute_high_order_solution_ss(A,b,viscL,mesh,...
            phys,quadrature,ev,dof_handler,opts,nonlin_opts,out_opts);
        
        % compute flux correction matrix
        F = flux_correction_matrix_ss(uH,DL-DH);
        
        % compute low-order solution
        uL = AL_mod \ b_mod;

        % prelimit flux corrections if user specified
        if (fct_opts.prelimit)
            F = prelimit_fluxes(F,uL);
        end
        
        % initialize solution iterate
        switch fct_opts.FCT_initialization
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
        for iter = 1:nonlin_opts.max_iter
            % compute limited flux correction sum
            [flim,Wminus,Wplus] = compute_limited_flux_sums_ss(uFCT,F,...
                AL_mod,b_mod,...
                sigma_min,sigma_max,source_min,source_max,mesh,phys,...
                dof_handler.n_dof,fct_opts,opts,...
                source_over_sigma_min, source_over_sigma_max);
            
            % plot
            if (out_opts.plot_FCT_iteration)
                plot_FCT(mesh.x,uH,uFCT,Wminus,Wplus,out_opts,true,true)
            end
            
            % compute correction rhs and modify for Dirichlet BC
            system_rhs = b + flim;
            if (opts.modify_for_strong_DirichletBC)
                system_rhs(1) = phys.inc;
            end
            
            % compute residual
            ss_res = system_rhs - AL_mod*uFCT;
            
            % test convergence of previous iteration
            converged = test_convergence(ss_res,nonlin_opts.nonlin_tol,iter);
            if (converged)
                fprintf('\tConverged at iteration %i\n',iter);
                break;
            end

            % solve modified system
            du = AL_mod \ ss_res;
            uFCT = uFCT + du;
        end
        
        if (~converged)
            error('FCT Solution did not converge');
        end
    else
        % compute initial conditions
        u_old = phys.IC(mesh.x);
        if phys.impose_BC_on_IC
          u_old(1) = phys.inc;
        end
        u_older = u_old;
        b_old = b;
        b_new = b_old;
        t = 0;
        time_step = 0;
        dt_old = dt_nominal;
        [DH,viscE] = compute_high_order_diffusion_matrix(u_old,...
            u_old,dt_old,viscL,mesh,phys,quadrature,ev,...
            dof_handler,opts.high_order_scheme);
        AH = A + DH;
        
        reached_end_of_transient = false;
        while (~reached_end_of_transient)
            % advance time step index
            time_step = time_step+1;
            
            % shorten time step if necessary
            if (t+dt_nominal >= opts.t_end)
                dt = opts.t_end - t;
                reached_end_of_transient = true;
            else
                dt = dt_nominal;
                reached_end_of_transient = false;
            end
            fprintf('Time step %i: t = %f->%f\n',time_step,t,t+dt);
            
            switch opts.temporal_scheme
                case 0 % steady-state
                    error('This should not be reachable');
                case 1 % explicit Euler
                    % perform high-order step
                    b = assemble_ss_rhs(t,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    [uH,DH] = high_order_step(u_older,u_old,dt_old,dt,...
                        A,b,MC,0,viscL,mesh,phys,quadrature,dof_handler,...
                        ev,opts);
                    
                    % perform FCT step
                    [uFCT,Wminus,Wplus] = FCT_step_explicit(u_old,uH,dt,...
                        ML,MC,AL,DH,DL,b,phys.inc,phys.speed,...
                        sigma_min,sigma_max,source_min,source_max,mesh,fct_opts,...
                        opts,phys.periodic_BC);
                case 2 % SSP3
                    % stage 1
                    u_old_stage = u_old;
                    u_older_stage = u_older;
                    t_stage = t;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,opts);
                    [uFCT_stage,Wminus,Wplus] = FCT_step_explicit(u_old_stage,uH_stage,dt,...
                        ML,MC,AL,DH,DL,b,phys.inc,phys.speed,...
                        sigma_min,sigma_max,source_min,source_max,mesh,fct_opts,...
                        opts,phys.periodic_BC);
                    
                   % stage 2
                    u_old_stage = uFCT_stage;
                    u_older_stage = u_older;
                    t_stage = t + dt;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,opts);
                    [uFCT_stage,Wminus,Wplus] = FCT_step_explicit(u_old_stage,uH_stage,dt,...
                        ML,MC,AL,DH,DL,b,phys.inc,phys.speed,...
                        sigma_min,sigma_max,source_min,source_max,mesh,fct_opts,...
                        opts,phys.periodic_BC);
                    
                    % stage 3
                    u_old_stage = 0.75*u_old + 0.25*uFCT_stage;
                    u_older_stage = u_older;
                    t_stage = t + 0.5*dt;
                    b = assemble_ss_rhs(t_stage,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    [uH_stage,DH] = high_order_step(u_older_stage,...
                        u_old_stage,dt_old,dt,A,b,MC,0,viscL,mesh,phys,...
                        quadrature,dof_handler,ev,opts);
                    [uFCT_stage,Wminus,Wplus] = FCT_step_explicit(u_old_stage,uH_stage,dt,...
                        ML,MC,AL,DH,DL,b,phys.inc,phys.speed,...
                        sigma_min,sigma_max,source_min,source_max,mesh,...
                        fct_opts,opts,phys.periodic_BC);
                    
                    % final combination
                    uFCT = 1/3*u_old + 2/3*uFCT_stage;
                case 3 % theta
                    % compute old steady-state residual
                    ss_res = b_old - AL*u_old;
                    
                    % compute system matrix
                    system_matrix = ML + opts.theta*dt*AL;
                    if (opts.modify_for_strong_DirichletBC)
                        system_matrix(1,:)=0; system_matrix(1,1)=1;
                    end
                    
                    % compute new ss rhs
                    b_new = assemble_ss_rhs(t+dt,quadrature,mesh,...
                        dof_handler,phys,opts.modify_for_weak_DirichletBC);
                    
                    % compute high-order solution
                    [uH,AH,DH] = compute_high_order_solution_theta(...
                        u_old,dt,MC,A,AH,b_old,b_new,viscL,quadrature,...
                        mesh,dof_handler,phys,ev,opts,...
                        nonlin_opts);

                    % compute flux correction matrix
                    F = flux_correction_matrix(u_old,uH,dt,DH,DL,MC,opts.theta);
                    
                    % compute theta ss rhs
                    b_theta = (1-opts.theta)*b_old + opts.theta*b_new;
                    
                    % compute solution
                    uL = compute_low_order_solution_theta(u_old,AL,...
                        ML,b_theta,dt,opts.theta,phys.inc,...
                        opts.modify_for_strong_DirichletBC);
                    
                    % initialize solution iterate
                    uFCT = zeros(dof_handler.n_dof,1);
                    
                    % iteration loop
                    for iter = 1:nonlin_opts.max_iter
                        % compute limited flux correction sum
                        [flim,Wminus,Wplus] = compute_limited_flux_sums(...
                            u_old,uFCT,uL,dt,...
                            ML,AL,b,F,sigma_min,sigma_max,source_min,...
                            source_max,mesh,opts.theta,dof_handler.n_dof,...
                            phys,fct_opts,opts);
                        
                        % compute system rhs
                        system_rhs = ML*u_old + (1-opts.theta)*dt*ss_res ...
                            + opts.theta*dt*b_new + dt*flim;
                        if (opts.modify_for_strong_DirichletBC)
                            system_rhs(1) = phys.inc;
                        end
                        
                        % compute residual
                        res = system_rhs - system_matrix*uFCT;
                        
                        % test convergence of previous iteration
                        converged = test_convergence(res,nonlin_opts.nonlin_tol,iter);
                        if (converged)
                            fprintf('\tConverged at iteration %i\n',iter);
                            break;
                        end
                        
                        % compute change in solution iterate
                        du = system_matrix \ system_rhs - uFCT;
    
                        % solve modified system
                        uFCT = uFCT + nonlin_opts.relax * du;

                        % plot
                        if (out_opts.plot_FCT_iteration)
                            figure;
                            clf;
                            hold on;
                            plot(mesh.x,uH,'rx');
                            plot(mesh.x,uL,'b+');
                            plot(mesh.x,uFCT,'g');
                            plot(mesh.x,Wminus,'k.');
                            plot(mesh.x,Wplus,'ko');
                            legend('High','Low','FCT(l+1)','W-(l)','W+(l)',...
                                'Location',out_opts.legend_location);
                            % pause
                            if (strcmp(out_opts.pause_type,'wait'))
                                waitforbuttonpress;
                            else
                                pause(out_opts.pausetime);
                            end
                        end
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
            reached_steady_state = ss_err < opts.ss_tol;
            if (reached_steady_state)
                fprintf('Transient terminated due to steady-state\n');
                break;
            end
            
            % plot
            if (out_opts.plot_FCT_transient)
                figure(1);
                clf;
                hold on;
                plot(mesh.x,uFCT,'g');
                plot(mesh.x,Wminus,'b--');
                plot(mesh.x,Wplus,'r--');
                legend('FCT','W-','W+','Location',out_opts.legend_location);
                % pause
                if (strcmp(out_opts.pause_type,'wait'))
                    waitforbuttonpress;
                else
                    pause(out_opts.pausetime);
                end
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
    u_exact = exact(xx,opts.t_end);
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
    plot(mesh.x,Wminus,'b--');
    plot(mesh.x,Wplus,'r--');
    legend_entries = char(legend_entries,'FCT','W-','W+');
end

% legend
legend(legend_entries,'Location',out_opts.legend_location);

%% Plot viscosity

% plot viscosity if requested
if (out_opts.plot_viscosity)
    figure; clf;
    
    % plot low-order viscosity if available
    if opts.low_order_scheme == 2 || opts.high_order_scheme ~= 1
        semilogy(mesh.x_center,viscL);
        legend_entries = char('Low-order viscosity');
    end
    
    hold on;
    
    % plot high-order viscosity if available
    if (opts.high_order_scheme ~= 1)
        semilogy(mesh.x_center,viscE,'x');
        legend_entries = char(legend_entries,'Entropy viscosity');
    end
    
    % legend
    legend(legend_entries,'Location',out_opts.legend_location);
end

%% Output

% save exact solution
if (save_exact_solution)
    exact_file = 'output/uexact.csv';
    dlmwrite(exact_file,[xx,u_exact],',');
end

% determine string to be appended to results for time discretization
switch opts.temporal_scheme
    case 0
        time_string = 'ss';
    case 1
        time_string = 'FE';
    case 2
        time_string = 'SSPRK33';
    case 3
        time_string = sprintf('theta%0.1f',opts.theta);
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
        low_order_file = ['output/uL_',time_string,'.csv'];
        dlmwrite(low_order_file,[mesh.x,uL_final],',');
    end
end

% determine string to be appended for high-order scheme
switch opts.high_order_scheme
    case 1
        high_order_string = 'Gal';
    case 2
        high_order_string = 'EV';
    case 3
        high_order_string = 'AltEV1';
    case 4
        high_order_string = 'AltEV2';
    otherwise
        error('Invalid high-order scheme chosen');
end

% save high-order solution
if (save_high_order_solution)
    if (compute_high_order)
        high_order_file = ['output/uH_',high_order_string,'_',time_string,'.csv'];
        dlmwrite(high_order_file,[mesh.x,uH_final],',');
    end
end

% save FCT solution
if (save_FCT_solution)
    if (compute_FCT)
        FCT_file = ['output/uFCT_',high_order_string,'_',time_string,'.csv'];
        dlmwrite(FCT_file,[mesh.x,uFCT],',');
    end
end

% save FCT bounds
if (save_FCT_bounds)
    if (compute_FCT)
        file = ['output/Umin_',high_order_string,'_',time_string,'.csv'];
        dlmwrite(file,[mesh.x,Wminus],',');
        file = ['output/Umax_',high_order_string,'_',time_string,'.csv'];
        dlmwrite(file,[mesh.x,Wplus],',');
    end
end

% save antidiffusion matrix
if (save_antidiffusion_matrix)
    if (compute_FCT)
        out_file = ['output/P.csv'];
        dlmwrite(out_file,F,',');
    end
end

% compute return values
return_value1 = 0;
return_value2 = 0;
switch return_value_option
    case 0 % nothing
        % do nothing
    case 1 % entropy residual
        % compute L^2 norm
        if (opts.is_steady_state)
            L2norm = evaluate_entropy_residual_l2norm_ss(...
                uH,mesh,phys,quadrature,ev,dof_handler);
        else
            error(['Entropy residual L^2 norm calculation only ',...
                'implemented for steady-state']);
        end
        
        % set return values
        return_value1 = L2norm;
        return_value2 = mesh.dx(1);
        
        % report log of mesh size and log of norm
        fprintf(['Note: main() function return values are:\n',...
            '  1. L2 norm of entropy residual\n',...
            '  2. mesh size\n']);
    case 2 % entropy jump
        % compute L^2 norm
        if (opts.is_steady_state)
            L2norm = evaluate_entropy_jump_l2norm_ss(...
                uH,mesh,phys,quadrature,ev,dof_handler);
        else
            error(['Entropy residual L^2 norm calculation only ',...
                'implemented for steady-state']);
        end
        
        % set return values
        return_value1 = L2norm;
        return_value2 = mesh.dx(1);
        
        % report log of mesh size and log of norm
        fprintf(['Note: main() function return values are:\n',...
            '  1. L2 norm of entropy residual\n',...
            '  2. mesh size\n']);
    otherwise
        error('Invalid return value option chosen');
end

end
