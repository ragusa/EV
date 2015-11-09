function euler_1d
clear;clc;clf;close all;

global hardcoded_disc
global data npar mesh

% OPTIONS
%================================================================================
hardcoded_disc = false;

% problem ID
% 1 = Leblanc's tube
% 2 = Sod's tube
data.pbID=2;

mesh.rand= false;
npar.nel = 50;
npar.porder = 1;               % FE polynomial degree
npar.qorder = npar.porder + 1; % quadrature order
% bypass my FEM disc and cell-dependent viscosity but use hardcoded linear
% with trapezoidal rule
if(hardcoded_disc)
    npar.porder=1; npar.qorder=2;
end
npar.temporal='ERK1';
npar.cfl =0.1;
npar.Cmax=0.5;%0.5;
npar.Cent=1.0;%1e0;
npar.lumpM = true;

time_approx_entropy = 'none';

plot_transient = false;
frequency=20; % how many time steps between transient plot outputs
%================================================================================

% number of degrees of freedom
npar.ndofs = npar.nel * npar.porder + 1 ;

% initialize data and parameters
init_data();
U_old = data.u0;
plot_vector(U_old,time,true)

% initialize old entropy values at previous time steps
S_m1 = zeros(npar.qorder,npar.nel);
S_m2 = S_m1;

% initialize time parameters
time = 0;
tend = data.tend;
dt   = npar.dt;
dt_old = -1;

iter=1;
tic
while time < tend

    if(time+dt>=tend)
        dt=tend-time;
        time = tend-dt+1e-15;
    end
    % time is the end time
    time=time+dt;
    fprintf('time beg=%9.5g, dt=%7.4d, time end=%9.5g \n',time-dt,dt,time);


    switch upper(npar.temporal)
        case{'ERK1'}
            % apply RK1 explicit
            k1 = ss_residual(U_old,time_approx_entropy);
            % value of the end of the time step
            U = U_old + dt*(npar.M\k1);
        case{'ERK2'}
            % apply RK2 explicit
            % first stage
            k1 = ss_residual(U_old,time_approx_entropy);
            U1 = U_old + (dt/2)*(npar.M\k1);
            % second stage
            k2 = ss_residual(U1,time_approx_entropy);
            % value of the end of the time step
            U = U_old + dt*(npar.M\k2);
        case{'ERK3'}
            % apply RK3 explicit
            % first stage
            k1 = ss_residual(U_old,time_approx_entropy);
            U1 = U_old + dt*(npar.M\k1);
            % second stage
            k2 = ss_residual(U1,time_approx_entropy);
            U2 = U_old + (dt/4)*(npar.M\(k1+k2));
            % third stage
            k3 = ss_residual(U2,time_approx_entropy);
            U = U_old + (dt/6)*(npar.M\(k1+k2+4*k3));
        case{'ERK4'}
            % apply RK4 explicit
            % first stage
            k1 = ss_residual(U_old,time_approx_entropy);
            U1 = U_old + (dt/2)*(npar.M\k1);
            % second stage
            k2 = ss_residual(U1,time_approx_entropy);
            U2 = U_old + (dt/2)*(npar.M\k2);
            % third stage
            k3 = ss_residual(U2,time_approx_entropy);
            U3 = U_old + dt*(npar.M\k3);
            % fourth stage
            k4 = ss_residual(U3,time_approx_entropy);
            U = U_old + (dt/6)*(npar.M\(k1+2*k2+2*k3+k4));
        otherwise
            error('unknown time integration %s',temporal);
    end

    % save old dt for bdf2
    npar.dt_old = dt;
    % get new dt
    dt = compute_new_dt(U);
    npar.dt = dt;
    % update how dsdt is computed in the entropy residual
    if(strcmpi(time_approx_entropy,'bdf1'))
        time_approx_entropy='bdf2';
    elseif(strcmpi(time_approx_entropy,'none'))
        time_approx_entropy='bdf1';
    end
    if (plot_transient)
        if(mod(iter,frequency)==0)
            plot_vector(U,time); drawnow
        end
    end
    % save old values for next time step
    U_old = U;
    iter=iter+1;
end

figure(2);
plot_vector(U,time,true);
toc

% % Uraw=U;
% % uu=smoothing(U,1);
% % 
% % plot_0=false;
% % if(plot_0)
% %     figure(1);
% %     plot(x_interp,uu_init,'r-', x,uu );
% % else
% %     figure(1);
% %     adr='R:\temporary\entropy\';
% %     adr='G:\RESEARCH-coupling\entropy\';
% %     plot(x_interp,uu,'.-b'); hold all;
% %     Z=load(sprintf('%sfort.10',adr)); plot(Z(:,1),Z(:,2),'.r-');
% %     plot(x_interp,uraw,'.-m');
% %     legend(['smo';'JLG';'raw']);
% %     hold off
% %     %     plot(x_interp,uu, x_interp,uraw,'r' );
% %     visc_=reshape(visc,qorder*nel,1);
% %     x_q=zeros(nel*qorder,1);
% %     for i=1:nel
% %         x1=x(i);  x2=x(i+1);
% %         aux=(x1+x2)/2 + (x2-x1)*xq/2;
% %         i1=(i-1)*qorder+1;
% %         i2=i1+qorder-1;
% %         x_q(i1:i2)=aux(1:end);
% %     end
% %     %     figure; plot(x_q,visc_);
% %     figure(2);
% %     subplot(3,1,1); plot(x_interp,uu);
% %     title(sprintf('Burgers 1d, porder=%i,tend=%g',porder,tend));
% %     subplot(3,1,2); plot(x_interp,uraw ); title('Raw solution -no averaging at tend')
% %     axis tight;
% %     subplot(3,1,3); plot(x_q,log10(visc_)); title('Viscosity at quadrature points')
% % 
% %     figure(3);
% %     plot(x_q,log10(visc_),'.-b'); hold all;
% %     Z=load(sprintf('%sfort.11',adr)); plot(Z(:,1),log10(Z(:,2)),'.r-');
% %     title('viscosity'); legend(['rag';'JLG']);
% %     hold off
% % 
% % end
% 
% % save sav_Burgers_RK3_99el_p1_CFL_0.5_partsF.mat x uu_init uu x_interp uraw x_q visc ;
% % save sav_Burgers_RK3_999el_p1_CFL_0.5_partsF.mat ;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original PDE:
%       du/dt + df/dx = 0
% Add stabilization through a viscous term:
%       du/dt + df/dx - d/dx(visc du/dx) = 0
%
% advection term
%   option a) int b df/dx = -int f db/dx
%   option b) df/dx =f' du/dx --> int b f' du/dx
% viscous term:
%   - int b d/dx(visc du/dx)= int visc du/dx db/dx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = ss_residual(U,temporal_approx)

global hardcoded_disc
global data npar mesh

% initialize ss residual vector
F=zeros(3*npar.ndofs,1);

% extract rho, mt, E at dofs
rho = U( 1:npar.ndofs               );
mt  = U((1:npar.ndofs) +  npar.ndofs);
E   = U((1:npar.ndofs) +2*npar.ndofs);
% compute pressure at dofs
P = data.eos(rho,mt,E);
% compute specific entropy at dofs
s = entropy(rho, P);
% compute velocity at dofs
v = mt./rho;
% compute speed of sound at dofs
c = sqrt(data.gamma*(P+data.Pinf)./rho);
% compute rieman invariants at dofsd
r = v+2*c/(data.gamma-1);
t = v-2*c/(data.gamma-1);

% normalization constants (use values at dofs, it is easier)
norm_c = max(abs(v+c - mean(v+c)));
norm_s = max(abs(s-mean(s)));

% compute the normalized jumps before computing the ss residual
% the jumps are calculated through a loop on all elements
[jumps,jump_vdsdx,jump_sdvdx] = compute_jumps(rho,mt,E,P,s,v,norm_s);

% rhs = -int_V { b div(f) )
%     =  int_V { grad(b).f } + BC terms (which we ignore b/c we have
%     Dirichlet BC for now)
% int_K grad(b).f = jac sum_q {w_q 1/jac dbdx_q f_q}
% for a given element, with linear FE and trapezoidal quadrature (w_q=1),
%    this becomes =  +/- 0.5 (fleft + fright)
%
% for iel=1:npar.nel
%     g=gn(iel,:);
%     % coef_u = uu(g(:)); % vector of length porder+1
%     % recall that shape is a matrix of size (qorder,porder+1)
%     % local_u = shape * coef_u; % vector of length qorder
%     % we do not need this for linear FE with tra
% end

if(hardcoded_disc)

    % compute speed
    v = mt./rho;
    v2= v.^2;
    % compute pressure
    P = data.eos(rho,mt,E);
    % compute speed of sound
    c = sqrt(data.gamma*(P+data.Pinf)./rho);
    % compute specific entropy
    s = entropy(rho, P);
    % specific internal energy
    e = E./rho - 0.5 * v2;
    % compute inviscid flux
    inviscid_flux = iflx(rho, mt, E, P, v);

    h=mesh.h;

    % cell gradients
    drhodx = diff(rho)/h;
    dvdx   = diff(v)/h;
    dedx   = diff(e)/h;
    dsdx   = diff(s)/h;
    % cell averages
    rhoa = rho(1:end-1)+0.5*diff(rho);
    va   = v(1:end-1)+0.5*diff(v);
    ea   = e(1:end-1)+0.5*diff(e);
    v2a  = v2(1:end-1)+0.5*diff(v2);
    rhova= mt(1:end-1)+0.5*diff(mt);

    % mass residual
    ind=1;
    for i=2:npar.ndofs-1
        F(i)= (inviscid_flux(i-1,ind)-inviscid_flux(i+1,ind))/2 + ...
            kappa(i)*drhodx(i) - kappa(i-1)*drhodx(i-1) ;
    end
    % mt residual
    ind=2;
    str=npar.stride(ind);
    for i=2:npar.ndofs-1
        F(i+str)= (inviscid_flux(i-1,ind)-inviscid_flux(i+1,ind))/2       + ...
            kappa(i)*va(i)  *drhodx(i) - kappa(i-1)*va(i-1)  *drhodx(i-1) + ...
            mu(i)   *rhoa(i)*dvdx(i)   - mu(i-1)   *rhoa(i-1)*dvdx(i-1)   ;
    end
    % energy residual
    ind=3;
    str=npar.stride(ind);
    for i=2:npar.ndofs-1
        F(i+str)= (inviscid_flux(i-1,ind)-inviscid_flux(i+1,ind))/2          + ...
            kappa(i)*rhoa(i) *dedx(i)   - kappa(i-1)*rhoa(i-1) *dedx(i-1)    + ...
            kappa(i)*ea(i)   *drhodx(i) - kappa(i-1)*ea(i-1)   *drhodx(i-1)  + ...
            -0.5*( kappa(i)*v2a(i)  *drhodx(i) - kappa(i-1)*v2a(i-1)  *drhodx(i-1) )+ ...
            mu(i)   *rhova(i)*dvdx(i)   - mu(i-1)   *rhova(i-1)*dvdx(i-1)    + ...
            kappa(i)*v2a(i)  *drhodx(i) - kappa(i-1)*v2a(i-1)  *drhodx(i-1)  ;
    end

else

    b  = npar.b; dbdx = npar.dbdx;
    dx = mesh.dx;
    wq = npar.wq;
    str= npar.stride(2);
    s_q_save_current = zeros(npar.qorder,npar.nel);
    s_q_m1 = zeros(npar.qorder,npar.nel);
    s_q_m2 = zeros(npar.qorder,npar.nel);
    r_q_m1 = zeros(npar.qorder,npar.nel);
    r_q_m2 = zeros(npar.qorder,npar.nel);
    t_q_m1 = zeros(npar.qorder,npar.nel);
    t_q_m2 = zeros(npar.qorder,npar.nel);

    for iel=1:npar.nel
        g=npar.gn(iel,:);
        dofs_rho = rho(g); % vector of length porder+1
        dofs_mt  = mt(g);  % vector of length porder+1
        dofs_E   = E(g);   % vector of length porder+1
        % not really dofs, but ...
        dofs_v = dofs_mt./dofs_rho;
        dofs_e = dofs_E./dofs_rho - 0.5 * dofs_v.^2;
        % recall that shape is a matrix of size (qorder,porder+1)
        rho_q = b * dofs_rho; % vector of length qorder
        mt_q  = b * dofs_mt;  % vector of length qorder
        E_q   = b * dofs_E;   % vector of length qorder
        e_q   = b * dofs_e;   % vector of length qorder
        % vector of length qorder
        drhodx_q = dbdx * dofs_rho;
        dvdx_q   = dbdx * dofs_v;
        dedx_q   = dbdx * dofs_e;
        % compute fluid speed
        v_q = mt_q./rho_q;
        % compute pressure
        P_q = data.eos(rho_q,mt_q,E_q);
        % compute speed of sound
        c_q = sqrt(data.gamma*abs(P_q+data.Pinf)./rho_q);

        % compute inviscid flux
        inviscid_flux_q = iflx(rho_q, mt_q, E_q, P_q, v_q);
        %         if(iel>28 & iel<31)
        %             rho_q
        %             mt_q
        %             E_q
        %             P_q
        %             v_q
        %             inviscid_flux_q
        %             disp('')
        %         end
        %-----------------------------------------------------------
        %%%%%%%%%%%%% entropy based viscous fluxes %%%%%%%%%%%%%%%%%
        %-----------------------------------------------------------
        % compute specific entropy
        s_q = entropy(rho_q, P_q);
        % alternate s_q
        s_q = b * s(g);
        s_q_save_current(:,iel) = s_q;
        % compute other riemannian invariants
        r_q = v_q + 2*c_q/(data.gamma-1); % wave at speed u-c
        t_q = v_q - 2*c_q/(data.gamma-1); % wave at speed u+c
        % alternate r_q
        r_q = b * r(g);
        t_q = b * t(g);
        r_q_save_current(:,iel) = r_q;
        t_q_save_current(:,iel) = t_q;
        % gradients of r,s,t
        dsdx_q = ( dbdx * s(g) ) * 2/mesh.dx(iel);
        drdx_q = ( dbdx * r(g) ) * 2/mesh.dx(iel);
        dtdx_q = ( dbdx * t(g) ) * 2/mesh.dx(iel);
        % time derivatives
        dsdt(:,iel) = visc_temporal_contribution(s_q,s_q_m1(:,iel),s_q_m2(:,iel),temporal_approx);
        drdt(:,iel) = visc_temporal_contribution(r_q,r_q_m1(:,iel),r_q_m2(:,iel),temporal_approx);
        dtdt(:,iel) = visc_temporal_contribution(t_q,t_q_m1(:,iel),t_q_m2(:,iel),temporal_approx);
        % residual for r,s,t
        Dsq = abs( dsdt(:,iel) + (v_q.*c_q).*dsdx_q ) / ( norm_s );
        Drq = abs( drdt(:,iel) + (v_q.*c_q).*drdx_q - c_q .* Dsq ) / norm_c;
        Dtq = abs( dtdt(:,iel) + (v_q.*c_q).*dtdx_q + c_q .* Dsq ) / norm_c;
        % entropy visc = Cent * h^2 * ( residual + jumps )
        mu_ent(:,iel) = ( npar.Cent * mesh.dx(iel)^2 ) * ...
            ( max(max(Dsq(:),Drq(:)),Dtq(:)) + max(jumps(:,iel)) );
        k_ent(:,iel) = mu_ent(:,iel);
        % first order visc = Cmax * h * speed
        mu_max = npar.Cmax * mesh.dx(iel) * max( v_q + c_q);
        k_max  = npar.Cmax * mesh.dx(iel) * max( v_q + c_q);
        % take minimum
        visc_ka = min (k_ent(:,iel),  k_max );
        visc_mu = min (mu_ent(:,iel), mu_max);
        % compute viscous flux
        viscous_flux_q = vflx(visc_ka, visc_mu, rho_q, drhodx_q, v_q, dvdx_q, e_q, dedx_q);
        %-----------------------------------------------------------
        %%% END END %%% entropy based viscous fluxes %%% END END %%%
        %-----------------------------------------------------------

        % compute total flux (add jacobian term to viscous flux)
        flux_q = inviscid_flux_q + 2/mesh.dx(iel) * viscous_flux_q;

        % caveat: dbdx(iq,itest_funct), so need to transpose here
        % to have a vector ordered along itest_func !
        % mass residual
        ind=1;
        F(g) = F(g) + dbdx' * (wq.*flux_q(:,ind)) ;
        % mt residual
        ind=2;
        g=g+str;
        F(g) = F(g) + dbdx' * (wq.*flux_q(:,ind)) ;
        % energy residual
        ind=3;
        g=g+str;
        F(g) = F(g) + dbdx' * (wq.*flux_q(:,ind)) ;

    end
    % save older entropies
    s_q_m2 = s_q_m1;
    s_q_m1 = s_q_save_current;
    r_q_m2 = r_q_m1;
    r_q_m1 = r_q_save_current;
    t_q_m2 = t_q_m1;
    t_q_m1 = t_q_save_current;
end

% apply Dirichlet to the ss residual
F(1          + npar.stride) = 0;
F(npar.ndofs + npar.stride) = 0;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pres=IGEOS(rho,mt,E)
% computes pressure for ideal gas law
% E is the total energy
global data

pres = (data.gamma-1)*( E - 0.5 * mt.^2 ./rho );

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pres=SGEOS(rho,mt,E)
% computes pressure for ideal gas law
% E is the total energy
global data

e = E./rho - 0.5 * ( mt./rho).^2; % specific internal energy
pres = (data.gamma-1)*( e-data.q ) .*rho  - data.gamma * data.Pinf;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = entropy(rho, P)

global data

% SPECIFIC entropy (``little'' s)
% The entropy is always to be computed with the real values of the density,
% momentum and total energy. Then in the function called pressure the area
% is set to 1.
s =  data.Cv * ( log(abs(P+data.Pinf)./rho.^(data.gamma)) ) + data.q_prime;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dt = compute_new_dt(U)
global mesh npar data

rho = U( 1:npar.ndofs              ) ;
mt  = U((1:npar.ndofs)+  npar.ndofs) ;
E   = U((1:npar.ndofs)+2*npar.ndofs) ;

% compute pressure
P = data.eos(rho,mt,E);
% compute speed of sound
csound = sqrt(data.gamma*abs(P+data.Pinf)./rho);
eig_max = max(csound+abs(mt./rho));

% compute the next time step
dt = mesh.h / (eig_max * npar.porder) * npar.cfl;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inviscid_flux = iflx(rho, mt, E, P, v )
% compute the inviscid flux
% E= total energy, v=speed
inviscid_flux = zeros(length(rho),3);

inviscid_flux(:,1) = mt;
inviscid_flux(:,2) = v.*mt+P;
inviscid_flux(:,3) = v.*(E+P);

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function viscous_flux = vflx(visc_ka, visc_mu, rho, drhodx, v, dvdx, e, dedx )
% compute the viscous flux
% e = specific internal energy, v=speed
viscous_flux = zeros(length(rho),3);

viscous_flux(:,1) = - visc_ka .* drhodx;
viscous_flux(:,2) = - visc_mu .* rho .* dvdx + v .* viscous_flux(:,1) ;
viscous_flux(:,3) = - visc_ka .* ( rho .* dedx + e .* drhodx ) ...
    -0.5 * v.^2 .* viscous_flux(:,1) + v .* viscous_flux(:,2) ;

return
end


%%%%%%%%%%%%%%%%%% JUMPS   %%%%%%%%%%%%%%%%%

function [jumps,varargout]=compute_jumps(rho,mt,E,P,s,v,norm_s)
% compute jumps (normal derivatives at interval extremities) of:
% d(u.s)/dx = u dsdx + s dudx

global hardcoded_disc
global data npar mesh

nel = npar.nel;

jumps = zeros(2,nel);

if(hardcoded_disc)

    dsdx  = zeros(nel,1);
    dudx  = zeros(nel,1);
    jump_vdsdx = zeros(2,nel);
    jump_sdvdx = zeros(2,nel);

    dsdx = diff(s)/mesh.h;
    dvdx = diff(v)/mesh.h;

    for i=2:nel-1
        jump_vdsdx(1,i) = v(i)  *( dsdx(i)-dsdx(i-1) );
        jump_vdsdx(2,i) = v(i+1)*( dsdx(i+1)-dsdx(i) );
        jump_sdvdx(1,i) = s(i)  *( dvdx(i)-dvdx(i-1) );
        jump_sdvdx(2,i) = s(i+1)*( dvdx(i+1)-dvdx(i) );
    end
    i=1;
    jump_vdsdx(2,i) = v(i+1)*( dsdx(i+1)-dsdx(i) );
    jump_sdvdx(2,i) = s(i+1)*( dvdx(i+1)-dvdx(i) );
    i=nel;
    jump_vdsdx(1,i) = v(i)  *( dsdx(i)-dsdx(i-1) );
    jump_sdvdx(1,i) = s(i)  *( dvdx(i)-dvdx(i-1) );


else

    dsdx_q  = zeros(2,nel);
    dudx_q  = zeros(2,nel);
    % s_q     = zeros(2,nel);
    % v_q     = zeros(2,nel);

    b    = npar.bJ;
    dbdx = npar.dbdxJ;

    v_K=zeros(nel+1,1);
    s_K=zeros(nel+1,1);
    for iel=1:nel,
        g=npar.gn(iel,:);
        % coef_u = uu(g(:)); % vector of length porder+1
        % recall that shapeJ is a matrix of size (2,porder+1)
        v_q = b * v(g); % vector of length 2
        s_q = b * s(g); % vector of length 2
        if(iel==1)
            v_K(1:2) = v_q;  s_K(1:2) = s_q;
        else
            v_K(iel+1) = v_q(2);
            s_K(iel+1) = s_q(2);
        end
        if(iel>1) %perform a check of continuity
            if( abs( v_K(iel)-v_q(1) ) >1e-8 ) , error('v_q not continuous???'); end
            if( abs( s_K(iel)-s_q(1) ) >1e-8 ) , error('s_q not continuous???'); end
        end
        dvdx_q(:,iel) = 2/mesh.dx(iel) * dbdx * v(g); % vector of length 2
        dsdx_q(:,iel) = 2/mesh.dx(iel) * dbdx * s(g); % vector of length 2
    end

    jump_vdsdx = aux_jump(v_K,dsdx_q,nel);
    jump_sdvdx = aux_jump(s_K,dvdx_q,nel);
    jumps = jump_vdsdx + jump_sdvdx;

end

% normalize jumps here
jump_vdsdx = jump_vdsdx / norm_s;
jump_sdvdx = jump_sdvdx / norm_s;
jumps = jumps / norm_s;

% extra output var
varargout{1} =  jump_vdsdx;
varargout{2} =  jump_sdvdx;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [jumps]=aux_jump(f_K,dgdx_q,nel,jumps)

if(nargin<4),
    jumps=zeros(2,nel);
end
for iel=2:nel-1
    jumps(:,iel) = jumps(:,iel) + ...
        max( abs( f_K(iel)  *(dgdx_q(1,iel)-dgdx_q(2,iel-1)) ),...
        abs( f_K(iel+1)*(dgdx_q(2,iel)-dgdx_q(1,iel+1)) ));
end
iel=1;
jumps(:,iel) = jumps(:,iel) + abs( f_K(iel+1)*(dgdx_q(2,iel)-dgdx_q(1,iel+1)) );
iel=nel;
jumps(:,iel) = jumps(:,iel) + abs( f_K(iel)  *(dgdx_q(1,iel)-dgdx_q(2,iel-1)) );

return
end

%%%%%%%%%%%%%% end of JUMPS   %%%%%%%%%%%%%%

function dsdt=visc_temporal_contribution(s_q, s_q_m1, s_q_m2, temporal_approx)

global npar

% temporal contribution to the entropy residual
if(strcmpi(temporal_approx,'none'))
    dsdt=zeros(length(s_q),1);
elseif(strcmpi(temporal_approx,'bdf1'))
    dsdt=( s_q - s_q_m1 ) /npar.dt;
else
    k0 = (2*npar.dt + npar.dt_old)/(npar.dt * (npar.dt + npar.dt_old));
    k1 = (npar.dt + npar.dt_old)/(npar.dt * npar.dt_old);
    k2 = npar.dt /(npar.dt_old * (npar.dt + npar.dt_old));
    dsdt = k0*s_q + k1*s_q_m1 + k2*s_q_m2;
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_vector(U,time,opts)
% Function to plot the density, the momentum, the pressure and the energy
% of the fluid
% init vectors

global mesh npar data

rho = U( 1:npar.ndofs              ) ;
mt  = U((1:npar.ndofs)+  npar.ndofs) ;
E   = U((1:npar.ndofs)+2*npar.ndofs) ;
% compute pressure
P = data.eos(rho,mt,E);
% compute speed of sound
csound = sqrt(data.gamma*abs(P+data.Pinf)./rho);
% fluid speed
v = mt./rho;
% mach number
Mach = abs(v./csound);
% temperature
T = abs(P+data.Pinf)./((data.gamma-1)*data.Cv*rho);
% entropy
s = entropy(rho, P);

% Space vector
X = mesh.x_interp;

% Plots
subplot(3,3,1);
plot(X,rho,'-'); title('Density'); xlabel('x(m)'); ylabel('Density(kg/m3)'); hold on;
%
subplot(3,3,2);
plot(X,mt,'-'); title('Momentum'); xlabel('x(m)'); ylabel('Momentum(kg/m2.s)'); hold on;
%
subplot(3,3,3);
plot(X,E,'-'); title('Energy'); xlabel('x(m)'); ylabel('Energy(J)'); hold on;
%
subplot(3,3,4);
plot(X,P,'-'); title('Pressure'); xlabel('x(m)'); ylabel('Pressure(Pa)'); hold on;
%
subplot(3,3,5);
plot(X,T,'-'); title('Temperature'); xlabel('x(m)'); ylabel('Temperature(K)'); hold on;
%
subplot(3,3,6);
plot(X,s,'-'); title('Entropy'); xlabel('x(m)'); ylabel('Entropy'); hold on;
%
subplot(3,3,7);
plot(X,v,'-'); title('Velocity'); xlabel('x(m)'); ylabel('Velocity(m/s)'); hold on;
%
subplot(3,3,8);
plot(X,Mach,'-'); title('Mach number');  xlabel('x(m)'); ylabel('Mach number'); hold on;
%
subplot(3,3,9);
plot(X,csound,'-'); title('Speed of sound'); xlabel('x(m)'); ylabel('Speed of sound'); hold on;

if(nargin<3 || opts==false), return; end

figure(3);

[ds,us,ps,es,xi]=ers(data.ql,data.qr,1,npar.nel,data.gamma,time);
% Plots
subplot(2,2,1);
plot(X,rho,'-.',xi,ds,'r-'); title('Density'); xlabel('x(m)'); ylabel('Density(kg/m3)'); hold on;
%
subplot(2,2,2);
plot(X,v,'-.',xi,us,'r-'); title('Velocity'); xlabel('x(m)'); ylabel('Velocity(m/s)'); hold on;
%
subplot(2,2,3);
plot(X,P,'-.',xi,ps,'r-'); title('Pressure'); xlabel('x(m)'); ylabel('Pressure(Pa)'); hold on;
%
subplot(2,2,4);
% ie=P./rho./(data.gamma-1);
plot(X,E,'-.',xi,ds.*es,'r-'); title('Energy'); xlabel('x(m)'); ylabel('Energy(J)'); hold on;
%

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function init_data()

global hardcoded_disc
global data npar

% select problem data
switch(data.pbID)
    case{1} % Leblanc shock tube
        data.L=1;
        data.membrane=data.L/2;
        % parameters for EOS and entropy
        data.eos = @IGEOS;
        data.gamma = 5/3; % gamma parameter
        data.Cv = 1 / (data.gamma-1); % heat capacity
        data.q       = 0; % q value
        data.Pinf    = 0; % Pinf
        data.q_prime = 20; % for positive entropy
        % boundary conditions and initial sound speeds
        v_left = 0;    v_rite = 0;
        rho_left = 1;  rho_rite = 0.001;
        mt_left = rho_left*v_left;
        mt_rite = rho_rite*v_rite;
        e_left = 0.1; e_rite = 1e-7;
        E_left = rho_left*e_left + 0.5*rho_left*v_left^2 ;
        E_rite = rho_rite*e_rite + 0.5*rho_rite*v_rite^2;
        %
        data.bc_value(:,1) = [rho_left; mt_left; E_left]; % x=0 (rho,m,E)
        data.bc_value(:,2) = [rho_rite; mt_rite; E_rite]; % x=L (rho,m,E)
        %data.bc_value(:,1) = [1;3;5]
        %data.bc_value(:,2) = [2;4;6]
        % initial speeds of sound
        Pl = data.eos(rho_left,mt_left,E_left);
        Pr = data.eos(rho_rite,mt_rite,E_rite);
        data.cl = sqrt( data.gamma * (Pl+data.Pinf) / rho_left);
        data.cr = sqrt( data.gamma * (Pr+data.Pinf) / rho_rite);
        data.ql=[1,    0,  0.1*2/3];
        data.qr=[0.001,0.0,0.0000001*2/3];

        data.tend=0.4;

    case{2} % Sod shock tube
        data.L=1;
        data.membrane=data.L/2;
        % parameters for EOS and entropy
        data.eos = @IGEOS;
        data.gamma = 1.4; % gamma parameter
        data.Cv = 1 / (data.gamma-1); % heat capacity
        data.q       = 0; % q value
        data.Pinf    = 0; % Pinf
        data.q_prime = 20; % for positive entropy
        % boundary conditions and initial sound speeds
        v_left = 0;    v_rite = 0;
        rho_left = 1;  rho_rite = 0.125;
        mt_left = rho_left*v_left;
        mt_rite = rho_rite*v_rite;
        p_left = 1.0;
        p_rite = 0.1;
        E_left = 2.5;
        E_rite = 0.25;
        %
        data.bc_value(:,1) = [rho_left; mt_left; E_left]; % x=0 (rho,m,E)
        data.bc_value(:,2) = [rho_rite; mt_rite; E_rite]; % x=L (rho,m,E)
        % initial speeds of sound
        Pl = data.eos(rho_left,mt_left,E_left);
        Pr = data.eos(rho_rite,mt_rite,E_rite);
        data.cl = sqrt( data.gamma * (Pl+data.Pinf) / rho_left);
        data.cr = sqrt( data.gamma * (Pr+data.Pinf) / rho_rite);
        data.ql=[rho_left, v_left, p_left];
        data.qr=[rho_rite, v_rite, p_rite];

        data.tend=0.25;

    otherwise
        error('unknown pbID %i',data.pbID);
end

% initialize spatial mesh
init_mesh();
%  initial condition
init_u0();

% store basis functions
[xq,wq] = GLNodeWt(npar.qorder);
[npar.b,npar.dbdx]=feshpln(xq,npar.porder);
% special shape functions for jumps
% h=0.01;
% [npar.bJ,npar.dbdxJ]=feshpln([-1 (-1+h) (1-h) 1],npar.porder);
[npar.bJ,npar.dbdxJ]=feshpln([-1 1],npar.porder);
if(hardcoded_disc)
    xq=[-1 1]; wq=[1 1];
    [npar.b ,npar.dbdx] =feshpln(xq,npar.porder);
    % [npar.bJ,npar.dbdxJ]=feshpln([-1 (-1+h) (1-h) 1],npar.porder);
    [npar.bJ,npar.dbdxJ]=feshpln([-1 1],npar.porder);
end
% store quadrature
npar.xq = xq;
npar.wq = wq;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function init_mesh()

global data mesh npar

% create the mesh
if(mesh.rand);
    npts = nel+1;
    if( ~isempty(data.membrane) ), npts = nel; end
    x = data.L * sort( rand(npts,1) );
    x = sort( [x data.membrane] );
    x(1) = 0; x(end) = data.L;
    mesh.x  = x'; clear x;
    mesh.dx = diff(mesh.x);
    mesh.h  = min(mesh.dx);
    error('disallowing random spatial mesh for now')
else
    x = linspace(0,data.L,npar.nel+1);
    if( ~isempty(data.membrane) )
        % a membrane is present
        % ind = find( abs(x-data.membrane)<1E-8 );
        % if(length(ind)==0), error('uniform mesh not aligned with membrane'); end
        ind = find( x<data.membrane );
        if(isempty(ind))
            error('membrane to the left of the domain');
        else
            i1=ind(end);
        end
        ind = find( x>data.membrane );
        if(isempty(ind))
            error('membrane to the right of the domain');
        else
            i2=ind(1);
        end
        switch (i2-i1)
            case{1}
                mesh.memb_elem = i1;
                x_middle=(x(i1)+x(i2))/2;
                if( abs(x_middle-data.membrane)<1E-8 )
                    [i1 x(i1) i2 x(i2) x_middle]
                    error('x_middle not aligned with membrane');
                end
            case{2}
                mesh.memb_elem = [i1 (i1+1)];
            otherwise
                error('i2(%i) should be i1(%i)+1 or 2',i2,i1);
        end
    end
    mesh.x  = x'; clear x;
    mesh.dx = diff(mesh.x);
    mesh.h  = min(mesh.dx);
end

% create a mesh of the nodal unknowns from the discretization
x_interp=zeros(npar.ndofs,1);
porder=npar.porder;
for i=1:npar.nel
    x1=mesh.x(i);
    x2=mesh.x(i+1);
    % this is not true with non-equally spaced Lagrangian poly for
    % example!!!!
    aux=linspace(x1,x2,porder+1);
    i1=(i-1)*porder+1;
    i2=i1+porder;
    x_interp(i1:i2-1)=aux(1:end-1);
end
x_interp(end)=mesh.x(end);
% store mesh for ploting
mesh.x_interp = x_interp; clear x_interp;

% first time step
safety=0.5;
npar.dt = safety * mesh.h * npar.cfl / max(data.cl,data.cr);

% connectivity
gn=zeros(npar.nel,(npar.porder+1));
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=[gn(iel-1,end) , gn(iel-1,2:end)+npar.porder ];
end
% if(logperiodic)
%     % periodic bc
%     gn(nel,end)=gn(1,1);
% end
npar.gn=gn;

% stride in between the variables rho, mt, E
stride=[0 npar.ndofs 2*npar.ndofs];
npar.stride=stride;

% compute local mass and stifness matrices
[npar.m,npar.k] = init_mat(npar.porder);

% create global mass matrix
nnz_one_var = ((porder+1)^2-1)*npar.nel+1;
M = spalloc(3*npar.ndofs,3*npar.ndofs,3*nnz_one_var);
for iel=1:npar.nel
    M(gn(iel,:),gn(iel,:)) = M(gn(iel,:),gn(iel,:)) + mesh.dx(iel)/2 * npar.m;
end
M(stride(2)+1:stride(2)+npar.ndofs,stride(2)+1:stride(2)+npar.ndofs) = M(1:npar.ndofs,1:npar.ndofs);
M(stride(3)+1:stride(3)+npar.ndofs,stride(3)+1:stride(3)+npar.ndofs) = M(1:npar.ndofs,1:npar.ndofs);
clear gn;
%
% apply Dirichlet BC
% -- left
ind=1;
M(ind+stride,:)=0;   M(:,ind+stride)=0;
for i=1:3
    M(ind+stride(i),ind+stride(i))=1;
end
% -- right
ind=npar.ndofs;
M(ind+stride,:)=0;   M(:,ind+stride)=0;
for i=1:3
    M(ind+stride(i),ind+stride(i))=1;
end
% keep M and inv of M
if(npar.lumpM)
    M=diag(sum(M));
end
npar.M  = M;
npar.iM = inv(M);
clear M;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function init_u0()
% compute the initial values

global data npar mesh

u0=zeros(3*npar.ndofs,1);

stride=npar.stride;
for iel=1:npar.nel
    gn=npar.gn(iel,:);
    len=length(gn);
    if( iel < mesh.memb_elem(1) )
        for i=1:len
            u0(gn(i)+stride)= data.bc_value(:,1);
        end
    elseif( iel > mesh.memb_elem(end) )
        for i=1:len
            u0(gn(i)+stride)= data.bc_value(:,2);
        end
    else
        if(length(mesh.memb_elem)==1)
            bc_slope=(data.bc_value(:,2)-data.bc_value(:,1))/(len-1);
            for i=1:len
                u0(gn(i)+stride)= (i-1)*bc_slope+data.bc_value(:,1);
                % hack
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                u0(gn(i)+stride)= data.bc_value(:,2);
            end
        else
            bc_middle= (data.bc_value(:,1)+data.bc_value(:,2))/2;
            if(iel==mesh.memb_elem(1))
                bc_slope=(bc_middle-data.bc_value(:,1))/(len-1);
                for i=1:len
                    u0(gn(i)+stride)= (i-1)*bc_slope+data.bc_value(:,1);
                    % hack
                    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if(i==len), u0(gn(i)+stride)= data.bc_value(:,2);end
                end
            else
                bc_slope2=(data.bc_value(:,2)-bc_middle)/(len-1);
                if(bc_slope~=bc_slope2), error('bc_slope/=bc_slope2'); end
                for i=2:len
                    u0(gn(i)+stride)= (i-1)*bc_slope2+bc_middle;
                end
            end
        end
    end
end
% figure(1);
% title('rho,mt,E @ initial conditions');
% subplot(3,1,1)
% plot(mesh.x_interp,u0(1:npar.ndofs),'.-')
% subplot(3,1,2)
% plot(mesh.x_interp,u0((1:npar.ndofs) +  npar.ndofs),'.-')
% subplot(3,1,3)
% plot(mesh.x_interp,u0((1:npar.ndofs) +2*npar.ndofs),'.-')

data.u0 = u0; clear u0;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = smoothing(uu,howmany)

out = uu;
for i=1:howmany
    out(1) = ( uu(end) + 4*uu(1) + uu(2) ) /6;
    out(2:end-1) = (uu(1:end-2) + 4*uu(2:end-1) + uu(3:end) ) /6;
    out(end) = (uu(end-1) + 4*uu(end) + uu(1) ) /6;
    uu=out;
end

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = smoothing_2(uu,howmany)

out = uu;
for i=1:howmany
    out(1) = max( [uu(end);uu(1:2)] ) ;
    for i=2:length(uu)-1
        out(i) = max( uu(i-1:i+1) );
    end
    out(end) = max( [uu(end-1:end);uu(1)] );
    uu=out;
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [m,k] = init_mat(p_order)

edof=p_order+1;
k=zeros(edof,edof);      % initialization of element matrix to zero
m=zeros(edof,edof);      % initialization of element matrix to zero
f=zeros(edof,1);         % initialization of rhs vector to zero

% numerical integration sampling points & weights
% mass matrix contains poly of order: pp = 2*p_order
% Gauss-Leg quad with q points can integrate poly of order 2*q+1
% thus, we need: 2*p_order = 2*q-1 => q=p_order+1
q=p_order+1;
[xq,w] = GLNodeWt(q);

% compute shape functions and derivatives at quad point
[shape,dhdr]=feshpln(xq,p_order);

for iq=1:q
    wtx=w(iq);
    %  compute element matrix
    for i=1:edof,
        for j=1:edof,
            k(i,j)=k(i,j)+( dhdr(iq,i)* dhdr(iq,j) )*wtx;
            m(i,j)=m(i,j)+(shape(iq,i)*shape(iq,j) )*wtx;
        end
        f(i)=f(i) + shape(iq,i)*wtx;
    end
end

impr=0;
if(impr>0)
    m
    k
    f
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [shapefun,dhdx]=feshpln (xv,p)

% EQUALLY-SPACED POLY
xd=linspace(-1,1,p+1);

shapefun=zeros(length(xv),p+1);
dhdx    =zeros(length(xv),p+1);

for i=1:p+1
    num=1.;
    den=1.;
    for j=1:p+1
        if(j~=i)
            num=num.*(xv-xd(j));
            den=den.*(xd(i)-xd(j));
        end
    end
    shapefun(:,i)=num./den;
end

for i=1:p+1
    sum=0.;
    den=1.;
    for j=1:p+1
        if(j~=i)
            num=1;
            for k=1:p+1
                if((k~=i)&&(k~=j))
                    num=num.*(xv-xd(k));
                end
            end
            sum=sum+num;
            den=den.*(xd(i)-xd(j));
        end
    end
    dhdx(:,i)=sum./den;
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = GLNodeWt(n)
% GLNodeWt  Nodes and weights for Gauss-Legendre quadrature of arbitrary order
%           obtained by solving an eigenvalue problem
%
% Synopsis:  [x,w] = GLNodeWt(n)
%
% Input:     n = order of quadrature rule
%
% Output:    x = vector of nodes
%            w = vector of weights

%  Algorithm based on ideas from Golub and Welsch, and Gautschi.  For a
%  condensed presentation see H.R. Schwarz, "Numerical Analysis: A
%  Comprehensive Introduction," 1989, Wiley.  Original MATLAB
%  implementation by H.W. Wilson and L.H. Turcotte, "Advanced Mathematics
%  and Mechanics Applications Using MATLAB," 2nd ed., 1998, CRC Press

beta   = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
J      = diag(beta,-1) + diag(beta,1);    % eig(J) needs J in full storage
[V,D]  = eig(J);
[x,ix] = sort(diag(D));  %  nodes are eigenvalues, which are on diagonal of D
w      = 2*V(1,ix)'.^2;  %  V(1,ix)' is column vector of first row of sorted V

return
end

%%%%%%%%%%%%%%%%%% JUMPS   %%%%%%%%%%%%%%%%%
% compute dF/dx and look at its jump
% F=entropy flux
% dF/dx = F' * du/dx = E' * f' * du/dx = f' * dE/dx
% use finite difference near extremities to get dE/dx
% since u is continuous, f'=df/du is continous at the interval extremities
% so, [[dF/dx]] = f'(\pm1) [[dE/dx]]

%%%%%%%%%%%%%% end of JUMPS   %%%%%%%%%%%%%%

% residual = dE/dt + dF/dx
% where F = entropy flux.
% however, dF/dx = F' * du/dx
% and F' = dF/du = E' * f'
% thus, residual = dE/dt + E'*f'*du/dx
%
% short demo: du/dt + df/dx = 0 <==> du/dt + f'(u) du/dx = 0
% for an entropy pair (E(u),F(u)) satisfying dE/dt + dF/dx = 0
% we have (smoothness condition) E'(u) du/dt + F'(u) du/dx = 0
% which is compared to: E'(u) [ du/dt + f'(u) du/dx ]
% % % % % % % %     figure(102);clf
% % % % % % % %     plot(resi(1,:)    ,'m.-'); hold all;
% % % % % % % %     plot(resi_JLG(1,:),'g.-'); legend(['rag';'jlg']); title('entro residual');  hold off
% % % % % % % %     clear resi;
% % % % % % % %     resi=resi_JLG(1,:)';
% % % % % % % %
% % % % % % % % if(c2>1.e5)
% % % % % % % %     visc = c1 * aux_dx_p .* speed;
% % % % % % % % else
% % % % % % % %     % here, the size h for the entropy visc is the minimum
% % % % % % % %     % between 2 quadrature points = dx/porder for lagrange polynomials
% % % % % % % %     % where the size for the maximum entropy (speed) is dx BUT
% % % % % % % %     % cmax(=c1) is a constant/porder in JLG paper, so it is like using
% % % % % % % %     % dx/porder everywhere in the formula
% % % % % % % %     visc1st= c1*aux_dx_p.*speed ;
% % % % % % % %     visc   = min( c1*aux_dx_p.*speed , c2*(aux_dx_quad).^2.*(resi+jumps) );
% % % % % % % %     visc_  = min( c1*aux_dx_p.*speed , c2*(aux_dx_quad).^2.*(jumps) );
% % % % % % % %     visc__ = min( c1*aux_dx_p.*speed , c2*(aux_dx_quad).^2.*(resi) );
% % % % % % % % end
% % % % % % % % figure(100); clf;
% % % % % % % % if(JLG_visc)
% % % % % % % %     if(c2<=1.e5)
% % % % % % % %         plot(log10(visc_ ),'g.-'); hold all;
% % % % % % % %         plot(log10(visc__),'m.-'); hold all;
% % % % % % % %     end
% % % % % % % %     plot(log10(visc  ),'b-'); hold all;
% % % % % % % %     plot(log10(visc1st ),'k.-'); legend(['J';'R';'B';'1']); hold off;
% % % % % % % % else
% % % % % % % %     if(c2<=1.e5)
% % % % % % % %         semilogy(max(visc1st),'k.-'); hold all;
% % % % % % % %         semilogy(max(visc_),'g.-'); hold all;
% % % % % % % %         semilogy(max(visc__),'m.-'); hold all;
% % % % % % % %     end
% % % % % % % % %     plot(log10(max(visc )),'b-'); legend(['J';'R';'B']); hold off;
% % % % % % % %     semilogy(max(visc ),'b-'); legend(['1';'J';'R';'B']); hold off;
% % % % % % % % end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact Riemann Solver for the 1D Euler equations,
% implemented from Toro's "Riemann Solvers and Numerical Methods for Fluid Dynamics"
% Ch. 4.9, pp. 153-162
% coded by Ivan Christov, Novemeber 2005, TAMU

%Leblanc test
% ql=[1,0,2/3] and qr=[0.001,0,0.0000001*2/3]
% timeout=6, domain =[0,9]

% functional implementation (the parameters are sent by the user)
% q = [rho,u,p]

function [ds,us,ps,es,xi] = ers(ql,qr,domlen,cells,gamma,timeout)
clear ds us ps es xi
mpa = 1;
diaph = domlen/2;
dl = ql(1); ul = ql(2); pl = ql(3);
dr = qr(1); ur = qr(2); pr = qr(3);

% gamma related constants
g1 = (gamma-1)/(2*gamma);
g2 = (gamma+1)/(2*gamma);
g3 = 2*gamma/(gamma-1);
g4 = 2/(gamma-1);
g5 = 2/(gamma+1);
g6 = (gamma-1)/(gamma+1);
g7 = (gamma-1)/2;
g8 = gamma-1;

% sound speeds
cl = sqrt(gamma*pl/dl);
cr = sqrt(gamma*pr/dr);

% check pressure's positivity
if g4*(cl+cr)<=(ur-ul)
    error('Initial data generates a vacuum!');
end

% get exact value of the pressure in the star region
% first, find educated guess of the value of the pressure in the star region
cup = 0.25*(dl+dr)*(cl+cr);
ppv = 0.5*(pl+pr)+0.5*(ul-ur)*cup;
ppv = max(0,ppv);
pmin = min(pl,pr);
pmax = max(pl,pr);
qmax = pmax/pmin;
quser = 2;
if (qmax<=quser) & ((pmin<=ppv) & (ppv<=pmax)) % select PVRS Riemann solver
    pm = ppv;
else
    if ppv<pmin % select two-rarefaction Riemann solver
        pq = (pl/pr)^g1;
        um = (pq*ul/cl+ur/cr+g4*(pq-1))/(pq/cl+1/cr);
        ptl = 1+g7*(ul-um)/cl;
        ptr = 1+g7*(um-ur)/cr;
        pm = 0.5*(pl*ptl^g3+pr*ptr^g3);
    else % select tow-shock Riemann solver with PVRS as estimate
        gel = sqrt((g5/dl)/(g6*pl+ppv));
        ger = sqrt((g5/dr)/(g6*pr+ppv));
        pm = (gel*pl+ger*pr-(ur-ul))/(gel+ger);
    end
end

% second, perform Newton iteration
niter = 5000;
tolpre = 10^(-12);
pold = pm;
for i=1:niter
    % evaluate the pressure function on the left and on the right
    if pold<=pl % rarefaction
        pratl = pold/pl;
        fl = g4*cl*(pratl^g1-1);
        fld = (1/(dl*cl))*pratl^(-g2);
    else % shock
        al = g5/dl;
        bl = g6*pl;
        qrtl = sqrt(al/(bl+pold));
        fl = (pold-pl)*qrtl;
        fld = (1-0.5*(pold-pl)/(bl+pold))*qrtl;
    end

    if pold<=pr % rarefaction
        pratr = pold/pr;
        fr = g4*cr*(pratr^g1-1);
        frd = (1/(dr*cr))*pratr^(-g2);
    else % shock
        ar = g5/dr;
        br = g6*pr;
        qrtr = sqrt(ar/(br+pold));
        fr = (pold-pr)*qrtr;
        frd = (1-0.5*(pold-pr)/(br+pold))*qrtr;
    end

    pm = pold-(fl+fr+ur-ul)/(fld+frd);
    change = 2*abs((pm-pold)/(pm+pold));
    if change<=tolpre break; end
    if pm<0 pm = tolpre; end
    pold = pm;
end

if i==niter
    warning('Newton-Raphson iteration did not achieve desired tolerance.');
end

% find the velocity in the star region
um = 0.5*(ul+ur+fr-fl);

% sample the solution
dx = domlen/cells;
for i=1:cells+1
    xi(i) = (i-1)*dx;%(i-0.5)*dx;
    if timeout==0
        if xi(i)<=diaph
            ds(i) = dl;
            us(i) = ul;
            ps(i) = pl;
        else
            ds(i) = dr;
            us(i) = ur;
            ps(i) = pr;
        end
    else
        s = (xi(i)-diaph)/timeout;
        if s<=um % left of contact discontinuity
            if pm<=pl % left rarefaction
                shl = ul-cl;
                if s<=shl % left data state
                    ds(i) = dl;
                    us(i) = ul;
                    ps(i) = pl;
                else
                    cml = cl*(pm/pl)^g1;
                    stl = um-cml;
                    if s>stl % star left state
                        ds(i) = dl*(pm/pl)^(1/gamma);
                        us(i) = um;
                        ps(i) = pm;
                    else % left fan
                        us(i) = g5*(cl+g7*ul+s);
                        cs = g5*(cl+g7*(ul-s));
                        ds(i) = dl*(cs/cl)^g4;
                        ps(i) = pl*(cs/cl)^g3;
                    end
                end
            else % left shock
                pml = pm/pl;
                sl = ul-cl*sqrt(g2*pml+g1);
                if s<=sl % left data state
                    ds(i) = dl;
                    us(i) = ul;
                    ps(i) = pl;
                else % left star state
                    ds(i) = dl*(pml+g6)/(pml*g6+1);
                    us(i) = um;
                    ps(i) = pm;
                end
            end
        else % right of contact discontinuity
            if pm>pr % right shock
                pmr = pm/pr;
                sr = ur+cr*sqrt(g2*pmr+g1);
                if s>=sr % right data state
                    ds(i) = dr;
                    us(i) = ur;
                    ps(i) = pr;
                else % star right state
                    ds(i) = dr*(pmr+g6)/(pmr*g6+1);
                    us(i) = um;
                    ps(i) = pm;
                end
            else % right rarefaction
                shr = ur+cr;
                if s>=shr % right data state
                    ds(i) = dr;
                    us(i) = ur;
                    ps(i) = pr;
                else
                    cmr = cr*(pm/pr)^g1;
                    str = um+cmr;
                    if s<=str % star right state
                        ds(i) = dr*(pm/pr)^(1/gamma);
                        us(i) = um;
                        ps(i) = pm;
                    else % right fan
                        us(i) = g5*(-cr+g7*ur+s);
                        cs = g5*(cr-g7*(ur-s));
                        ds(i) = dr*(cs/cr)^g4;
                        ps(i) = pr*(cs/cr)^g3;
                    end
                end
            end
        end
    end
end
es = ps./g8+ds.*(us.^2)./2;
% DEBUG INFORMATION
subplot(2,2,1); plot(xi,ds); ylabel('\rho'); xlabel('x');
subplot(2,2,2); plot(xi,us); ylabel('u'); xlabel('x');
subplot(2,2,3); plot(xi,ps); ylabel('p'); xlabel('x');
subplot(2,2,4); plot(xi,ps./ds./g8); ylabel('e'); xlabel('x');
return
end
