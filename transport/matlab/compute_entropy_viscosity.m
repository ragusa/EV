function viscE = compute_entropy_viscosity(...
    u_old,u_new,x,mu,sigma,source,dt,...
    v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,periodic_BC)

n = length(u_old); % number of dofs

if periodic_BC
    nel = n;  % number of elements
else
    nel = n-1;
end

g = [linspace(1,nel,nel)' linspace(2,nel+1,nel)']; % connectivity array
if periodic_BC
    g(end,:) = [ g(end,1) 1 ];
end

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
n_faces=n;
jump = zeros(n_faces,1);

for f = 1:n_faces % loop over interior faces
    if ~periodic_BC
        if f==1 || f==n_faces
            continue
        else % interior face when Dirchlet
            iL = f-1;   % left  cell index
            iR = f; % right cell index
        end
    else
        if f==1 
            iL = nel;   % left  cell index
            iR = f; % right cell index
        else % interior face when Dirchlet
            iL = f-1;   % left  cell index
            iR = f; % right cell index
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
    R = (E_new-E_old)/dt + mu*dEdx_new +...
        dEdu_new.*(sigma_q.*u_new_local - source_q);
    R_max = max(abs(R));
    
    % compute max jump
    iel_next = iel+1;
    if periodic_BC && iel==nel
        iel_next=1;
    end
    jump_max = max(jump(iel),jump(iel_next));
    
    % compute entropy viscosity
    viscE(iel) = (cE*R_max + cJ*jump_max)/ E_dev_max;
end