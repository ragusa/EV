function viscE = compute_entropy_viscosity(...
    u_old,u_new,dt,mesh,phys,quadrature,ev,dof_handler)

% unpack quadrature
zq   = quadrature.zq;
wq   = quadrature.wq;
v    = quadrature.v;
dvdz = quadrature.dvdz;
Jac  = quadrature.Jac;

% unpack physics
mu     = phys.mu;
sigma  = phys.sigma;
source = phys.source;
speed  = phys.speed;
periodic_BC = phys.periodic_BC;

% unpack mesh
x   = mesh.x;
nel = mesh.n_cell;

% unpack dof handler
n_dof = dof_handler.n_dof;
g     = dof_handler.connectivity;

% unpack entropy viscosity
cE            = ev.cE;
cJ            = ev.cJ;
entropy       = ev.entropy;
entropy_deriv = ev.entropy_deriv;


% compute domain average of entropy
E_integral = 0;
for iel = 1:nel
    u_new_local = v' * u_new(g(iel,:)); % solution evaluated at quad. points
    Eq = entropy(u_new_local);
    E_integral = E_integral + dot(wq,Eq)*Jac(iel);
end
L = x(end)-x(1);        % length of domain
E_avg = E_integral / L; % domain average of entropy

% compute max entropy deviation in domain for normalization constant
E_dev_max = 0;
for iel = 1:nel
    u_new_local = v' * u_new(g(iel,:)); % solution evaluated at quad. points
    Eq = entropy(u_new_local);
    E_dev_max = max(E_dev_max, max(abs(Eq-E_avg)));
end

% compute entropy jumps
n_face = n_dof;
jump = zeros(n_face,1);

for f = 1:n_face % loop over interior faces
    if periodic_BC
        if f == 1
            iL = nel; % left  cell index
            iR = f;   % right cell index
        else
            iL = f-1; % left  cell index
            iR = f;   % right cell index
        end
    else
        if f == 1 || f == n_face
            % jumps on exterior faces considered zero
            continue
        else
            iL = f-1; % left  cell index
            iR = f;   % right cell index
        end
    end
    uF = u_new(g(iL,2)); % solution on face
    dudx_L = dvdz' * u_new(g(iL,:)) / Jac(iL);
    dudx_R = dvdz' * u_new(g(iR,:)) / Jac(iR);
    dEdx_L = entropy_deriv(uF)*dudx_L(end);
    dEdx_R = entropy_deriv(uF)*dudx_R(1);
    jump(f) = abs(mu*(dEdx_L - dEdx_R));
end

% compute entropy viscosity
viscE = zeros(nel,1);
for iel = 1:nel
    % evaluate new and old solution at quadrature points
    u_new_local_nodes = u_new(g(iel,:));
    u_old_local_nodes = u_old(g(iel,:));
    u_new_local = v' * u_new_local_nodes;
    u_old_local = v' * u_old_local_nodes;
    dudx_new_local  = dvdz' * u_new(g(iel,:)) / Jac(iel);
    % evaluate entropy and derivative
    E_new = entropy(u_new_local);
    E_old = entropy(u_old_local);
    dEdu_new = entropy_deriv(u_new_local);
    dEdx_new = dEdu_new.*dudx_new_local;
    % compute quadrature point positions
    xq = get_quadrature_point_positions(x,iel,zq);
    % evaluate cross section and source at quadrature points
    sigma_q  = sigma(xq);
    source_q = source(xq);
    
    % compute maximum entropy residual at quadrature points
    R = (E_new-E_old)/(speed*dt) + mu*dEdx_new +...
        dEdu_new.*(sigma_q.*u_new_local - source_q);
    R_max = max(max(0,R));
    
    % compute max jump
    faceL = g(iel,1);
    faceR = g(iel,2);
    jump_max = max(jump(faceL),jump(faceR));
    
    % compute entropy viscosity
    if abs(E_dev_max) < 1.0e-100
       viscE(iel) = 0.0;
    else
       viscE(iel) = (cE*R_max + cJ*jump_max) / E_dev_max;
    end
end

% smooth the entropy viscosity if specified
if ev.smooth_entropy_viscosity
    % get the smoothing weight
    weight = ev.smoothing_weight;

    % save the unsmoothed entropy viscosity
    viscE_uns = viscE;

    % first cell
    viscE(1) = (viscE_uns(1)*weight + viscE_uns(2))/(1+weight);

    % interior cells
    for iel = 2:nel-1
        viscE(iel) = (viscE_uns(iel)*weight + viscE_uns(iel-1) + ...
            viscE_uns(iel+1)) / (2+weight);
    end

    % first cell
    viscE(nel) = (viscE_uns(nel)*weight + viscE_uns(nel-1))/(1+weight);
end
