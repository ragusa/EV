close all; clc;

% ML*du/dt-(K+D)*u = b
% MC*du/dt-(K  )*u = b

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
impose_weakly = 1; % 0 = impose Dirichlet strongly
                   % 1 = impose Dirichlet weakly
%==========================================================================

n_dof = nel + 1; % number of dofs

x=linspace(0,len,nel+1); % x points for plotting

% build matrices
[MC,K,ML,D,b]=build_matrices(len,nel,omega,sigma,src,inc,impose_weakly);
ML = diag(ML);
B = (ML-MC)/ML; % B = (ML-MC)*ML^-1 (Guermond)
AL = -(K+D);    % low-order steady-state matrix
AH = -K;        % high-order steady-state matrix
% compute modified lumped mass matrix for Dirichlet BC
ML_mod = ML;
ML_mod(1,:)=0;ML_mod(1,1)=1;
% compute modified consistent mass matrix for Dirichlet BC
MC_mod = MC;
MC_mod(1,:)=0;MC_mod(1,1)=1;

% compute dt using CFL condition
dt = CFL*ML(1,1)/AL(1,1);
for i = 2:n_dof
    dt = min(dt, CFL*ML(i,i)/AL(i,i));
end

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
    rhs = ML*u_old + dt*(b-AL*u_old);
    if (impose_weakly) % impose the Dirichlet BC weakly
        u_new = ML \ rhs;
    else % impose the Dirichlet BC strongly
        % modify rhs for Dirichlet BC
        rhs(1) = inc;
        % solve modified system
        u_new = ML_mod \ rhs;
    end
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
    rhs = MC*u_old + dt*(b-AH*u_old);
    if (impose_weakly) % impose the Dirichlet BC weakly
        u_new = MC \ rhs;
    else % impose the Dirichlet BC strongly
        % modify rhs for Dirichlet BC
        rhs(1) = inc;
        % solve modified system
        u_new = MC_mod \ rhs;
    end
    % reset u_old
    u_old = u_new;
end
plot(x,u_new,'b-+');
legend_entries = char(legend_entries,'High-order');

% FCT solve, incorrect correction
%==========================================================================
fprintf('Computing FCT (incorrect correction) solution...\n');
u_old = u0;
for n = 1:n_dt
    % low-order solve
    %----------------
    % compute rhs
    rhs = ML*u_old + dt*(b-AL*u_old);
    if (impose_weakly) % impose the Dirichlet BC weakly
        uL = ML \ rhs;
    else % impose the Dirichlet BC strongly
        % modify rhs for Dirichlet BC
        rhs(1) = inc;
        % solve modified system
        uL = ML_mod \ rhs;
    end

    % FCT solve
    %----------------
    % compute flux correction matrix. You can check the following
    %    sum(F,2) = -dt*D*u_old + dt*B*(K*u_old + b), where sum(F,2) is row-sum
    correction = compute_flux_correction_matrix(u_old,dt,D,B,K,b);
    % compute the upper and lower bound for max principle
    [W_max,W_min] = compute_max_principle_bounds(u_old,dt,ML,AL,b);
    % compute limiting coefficients
    switch limiting_option
        case 0 % no correction
            limiter = zeros(n_dof,n_dof);
        case 1 % full correction (no limiting)
            limiter = ones(n_dof,n_dof);
        case 2 % normal limiting
            limiter = compute_limiting_coefficients(correction,uL,ML,W_max,W_min);
        otherwise
            error('Invalid limiting option');
    end
    % compute correction rhs
    rhs = ML*uL + sum((limiter.*correction),2);
    if (impose_weakly) % impose the Dirichlet BC weakly
        u_new = ML \ rhs;
    else % impose the Dirichlet BC strongly
        % modify rhs for Dirichlet BC
        rhs(1) = inc;
        % solve modified system
        u_new = ML_mod \ rhs;
    end
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
if (impose_weakly)
    BC_error_string = '';
else 
    BC_error_string = ', incorrect';
end
FCT_legend_string = ['FCT, ',limiter_string,BC_error_string];
legend_entries = char(legend_entries,FCT_legend_string);
plot(x,W_min,'m:');
legend_entries = char(legend_entries,'FCT Min bound');
plot(x,W_max,'m:');
legend_entries = char(legend_entries,'FCT Max bound');

if (~impose_weakly)
    % FCT solve, correct correction
    %==========================================================================
    fprintf('Computing FCT (correct correction) solution...\n');
    u_old = u0;
    for n = 1:n_dt
        % low-order solve
        % compute rhs
        rhs = ML*u_old + dt*(b-AL*u_old);
        % modify rhs for Dirichlet BC
        rhs(1) = inc;
        % solve modified system
        uL = ML_mod \ rhs;
        
        % consult my notes FCT.pdf, Equation (12):
        %    correction_rowsum = \mathcal{A}*1
        Ftr = ML*u_old + dt*(b-AL*u_old);
        Ftr(1) = inc;
        Gtr = MC*u_old + dt*(b-AH*u_old);
        Gtr(1) = inc;
        correction_rowsum = ML/MC_mod*Gtr - ML/ML_mod*Ftr;
        u_new = uL + ML\correction_rowsum;
        
        % reset u_old
        u_old = u_new;
    end
    plot(x,u_new,'m-o');
    legend_entries = char(legend_entries,'FCT, correct');
end

legend(legend_entries);