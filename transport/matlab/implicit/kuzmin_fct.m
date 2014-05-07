close all; clc;
% Solves a transient transport equation:
%   MC*du/dt-(K  )*u = b (high-order system)
%   ML*du/dt-(K+D)*u = b (low-order system)
% using either explicit Euler or implicit Euler

len=10;    % length of domain
nel=10;    % number of elements
omega=1.0; % omega
sigma=1;   % sigma
src=0;     % source
inc=1;     % incoming flux
n_dt  = 10;      % number of time steps
CFL = 0.8;       % CFL number
limiting_option = 2; % 0 = set limiting coefficients to 0 (no correction)
                     % 1 = set limiting coefficients to 1 (full correction)
                     % 2 = compute limiting coefficients normally
use_explicit = 0; % 0 = use implicit Euler
                  % 1 = use explicit Euler
%==========================================================================

n_dof = nel + 1; % number of dofs

x=linspace(0,len,nel+1); % x points for plotting

% build matrices
[MC,ML,K,D,b]=build_matrices(len,nel,omega,sigma,src);
ML = diag(ML);
A = -(K+D);    % low-order steady-state system matrix

% compute dt using CFL condition, even for implicit, just so that time
% steps will be equal between explicit and implicit, for comparison
dt = CFL*ML(1,1)/A(1,1);
for i = 2:n_dof
    dt = min(dt, CFL*ML(i,i)/A(i,i));
end

% define low-order and high-order transient sytem matrices to be inverted
if (use_explicit)
    AL = ML;
    AH = MC;
else
    AL = ML-dt*(K+D);
    AH = MC-dt*K;
end
% compute modified transient sytem matrix for Dirichlet BC
AL_mod = AL;
AL_mod(1,:)=0; AL_mod(1,1)=1;
% compute modified transient sytem matrix for Dirichlet BC
AH_mod = AH;
AH_mod(1,:)=0; AH_mod(1,1)=1;

% initial conditions for pseudotransient (equal to exact steady-state solution)
if(sigma>eps)
    u0=src/sigma+(inc-src/sigma)*exp(-sigma*x'/omega);
else
    u0=inc+src/omega*x';
end

hold all;

% exact solution
%==========================================================================
xx=linspace(0,len,1000);
if(sigma>eps)
    exact=src/sigma+(inc-src/sigma)*exp(-sigma*xx/omega);
else
    exact=inc+src/omega*xx;
end
plot(xx,exact,'k-');
legend_entries = char('Exact');

% low-order solve
%==========================================================================
fprintf('Computing low-order solution...\n');
u_old = u0;
for n = 1:n_dt
    % compute rhs
    if (use_explicit)
        rhs = (ML - dt*A)*u_old + dt*b;
    else
        rhs = ML*u_old + dt*b;
    end
    % modify rhs for Dirichlet BC
    rhs(1) = inc;
    % solve modified system
    u_new = AL_mod \ rhs;
    % reset u_old
    u_old = u_new;
end
plot(x,u_new,'r-s');
legend_entries = char(legend_entries,'Low-order');

% high-order solve
%==========================================================================
fprintf('Computing high-order solution...\n');
u_old = u0;
for n = 1:n_dt
    % compute rhs
    if (use_explicit)
        rhs = (MC + dt*K)*u_old + dt*b;
    else
        rhs = MC*u_old + dt*b;
    end
    % modify rhs for Dirichlet BC
    rhs(1) = inc;
    % solve modified system
    u_new = AH_mod \ rhs;
    % reset u_old
    u_old = u_new;
end
plot(x,u_new,'b-+');
legend_entries = char(legend_entries,'High-order');

% FCT solve
%==========================================================================
fprintf('Computing FCT solution...\n');
u_old = u0;
for n = 1:n_dt
    % high-order solve
    %----------------
    % compute rhs
    if (use_explicit)
        rhs = (MC + dt*K)*u_old + dt*b;
    else
        rhs = MC*u_old + dt*b;
    end
    % modify rhs for Dirichlet BC
    rhs(1) = inc;
    % solve modified system
    uH = AH_mod \ rhs;

    % FCT solve
    %----------------
    % compute flux correction matrix
    F = compute_kuzmin_flux_correction_matrix(u_old,uH,dt,D,MC,use_explicit);
    % compute the upper and lower bound for max principle
    [W_max,W_min] = compute_max_principle_bounds(u_old,dt,ML,A,b);
    % compute limiting coefficients
    switch limiting_option
        case 0 % no correction
            limiter = zeros(n_dof,n_dof);
        case 1 % full correction (no limiting)
            limiter = ones(n_dof,n_dof);
        case 2 % normal limiting
            limiter = compute_limiting_coefficients_kuzmin(F,u_old,ML,W_max,W_min);
        otherwise
            error('Invalid limiting option');
    end
    % compute correction rhs
    if (use_explicit)
        rhs = (ML - dt*A)*u_old + dt*b + sum((limiter.*F),2);
    else
        rhs = ML*u_old + dt*b + sum((limiter.*F),2);
    end
    % modify rhs for Dirichlet BC
    rhs(1) = inc;
    % solve modified system
    u_new = AL_mod \ rhs;
    % check that max principle is satisfied
    check_max_principle(u_new,W_max,W_min);
    % reset u_old
    u_old = u_new;
end
plot(x,u_new,'g-x');
% create legend string for this set
switch limiting_option
    case 0 % no correction
        limiter_string = 'no correction';
    case 1 % full correction (no limiting)
        limiter_string = 'not limited';
    case 2 % normal limiting
        limiter_string = 'limited';
    otherwise
        error('Invalid limiting option');
end
FCT_legend_string = ['FCT, ',limiter_string];
legend_entries = char(legend_entries,FCT_legend_string);
plot(x,W_min,'m:');
legend_entries = char(legend_entries,'FCT Min bound');
plot(x,W_max,'m:');
legend_entries = char(legend_entries,'FCT Max bound');

legend(legend_entries);