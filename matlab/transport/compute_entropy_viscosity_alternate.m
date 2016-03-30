function viscE = compute_entropy_viscosity_alternate(...
    u_old,u_new,dt,mesh,phys,quadrature,ev,dof_handler,high_order_scheme)

if (high_order_scheme == 2)
    error('This function is not for normal entropy viscosity');
end

% unpack quadrature
nq   = quadrature.nq;
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

% option to use local entropy normalization
use_local_entropy_normalization = true;

% compute domain average of entropy
E_integral = 0;
for iel = 1:nel
    % evaluate solution at quadrature points
    u_new_local = v' * u_new(g(iel,:));

    % make transformation: y(u) = u*exp(sigma*x)
    xq = get_quadrature_point_positions(x,iel,zq);
    sigma_q  = sigma(xq);

    % compute entropy and add to integral
    y_new_local = u_new_local.*exp(xq.*sigma_q);
    Eq = entropy(y_new_local);
    E_integral = E_integral + dot(wq,Eq)*Jac(iel);
end
L = x(end)-x(1);        % length of domain
E_avg = E_integral / L; % domain average of entropy

% compute max entropy deviation in domain for normalization constant
% E_dev_max = 0;
% for iel = 1:nel
%     % evaluate solution at quadrature points
%     u_new_local = v' * u_new(g(iel,:));
% 
%     % make transformation: y(u) = u*exp(sigma*x)
%     xq = get_quadrature_point_positions(x,iel,zq);
%     sigma_q  = sigma(xq);
%     y_new_local = u_new_local.*exp(xq.*sigma_q);
% 
%     % compute entropy and deviation from average
%     if (high_order_scheme == 3)
%         cnorm = entropy(y_new_local) - E_avg;
%     else
%         cnorm = 0.5 * u_new_local .* u_new_local - exp(-2.0*sigma_q.*xq)*E_avg;
%     end
%      E_dev_max = max(E_dev_max, max(abs(cnorm)));
% end

% % compute entropy jumps
% n_face = n_dof;
% jump = zeros(n_face,1);
% 
% for f = 1:n_face % loop over interior faces
%     if periodic_BC
%         if f == 1
%             iL = nel; % left  cell index
%             iR = f;   % right cell index
%         else
%             iL = f-1; % left  cell index
%             iR = f;   % right cell index
%         end
%     else
%         if f == 1 || f == n_face
%             % jumps on exterior faces considered zero
%             continue
%         else
%             iL = f-1; % left  cell index
%             iR = f;   % right cell index
%         end
%     end
%     % get solution on face
%     uF = u_new(g(iL,2));
% 
%     % compute solution and gradients at quadrature points
%     u_new_L = v' * u_new(g(iL,:));
%     u_new_R = v' * u_new(g(iR,:));
%     dudx_L = dvdz' * u_new(g(iL,:)) / Jac(iL);
%     dudx_R = dvdz' * u_new(g(iR,:)) / Jac(iR);
% 
%     % make transformation: y(u) = u*exp(sigma*x)
%     sigma_F = sigma(x(f)); % compute cross section on face
%     yF = uF*exp(x(f)*sigma_F);
% 
%     % make transformation: y(u) = u*exp(sigma*x)
%     xq_L = get_quadrature_point_positions(x,iL,zq);
%     xq_R = get_quadrature_point_positions(x,iR,zq);
%     sigma_q_L  = sigma(xq_L);
%     sigma_q_R  = sigma(xq_R);
%     y_new_L = u_new_L.*exp(sigma_q_L.*xq_L);
%     y_new_R = u_new_R.*exp(sigma_q_R.*xq_R);
%     dydx_L = dudx_L .* exp(sigma_q_L .* xq_L) + ...
%       sigma_q_L .* exp(sigma_q_L .* xq_L) .* u_new_L;
%     dydx_R = dudx_R .* exp(sigma_q_R .* xq_R) + ...
%       sigma_q_R .* exp(sigma_q_R .* xq_R) .* u_new_R;
% 
%     % compute jumps
%     dEdx_L = entropy_deriv(yF)*dydx_L(end);
%     dEdx_R = entropy_deriv(yF)*dydx_R(1);
%     jump(f) = abs(mu*(dEdx_L - dEdx_R));
% end

% compute entropy viscosity
viscE = zeros(nel,1);
for iel = 1:nel
    % get quadrature point positions and compute cross section and source
    xq = get_quadrature_point_positions(x,iel,zq);
    sigma_q  = sigma(xq);

    % evaluate new and old solution at quadrature points
    u_new_local_nodes = u_new(g(iel,:));
    u_new_local = v' * u_new_local_nodes;
    dudx_new_local  = dvdz' * u_new(g(iel,:)) / Jac(iel);

    % computeentropy residual at quadrature points
    if (high_order_scheme == 3)
        % make transformation: y(u) = u*exp(sigma*x)
        y_new_local = u_new_local.*exp(sigma_q.*xq);
        dydx_new_local = dudx_new_local .* exp(sigma_q .* xq) + ...
            sigma_q .* exp(sigma_q .* xq) .* u_new_local;
        dEdy_new = entropy_deriv(y_new_local);
        
        R = dEdy_new.*dydx_new_local;
    elseif (high_order_scheme == 4)
        R = u_new_local .* dudx_new_local + sigma_q.*u_new_local.*u_new_local;
    end
    
%     % compute max jump
%     faceL = g(iel,1);
%     faceR = g(iel,2);
%     jump_max = max(jump(faceL),jump(faceR));
    
    % compute entropy viscosity
    if (use_local_entropy_normalization)
        if (high_order_scheme == 3)
            % compute normalization at each quadrature point
            normalization = entropy(y_new_local) - E_avg;
        elseif (high_order_scheme == 4)
            % compute normalization at each quadrature point
            normalization = 0.5 * u_new_local .* u_new_local ...
                - exp(-2*sigma_q.*xq)*E_avg;
        end

        % compute max entropy viscosity at quadrature points on cell
        viscE(iel) = 0.0;
        for q = 1:nq
            if abs(normalization(q)) < 1.0e-100
               viscEq = 0.0;
            else
               viscEq = abs(cE*R(q) / normalization(q));
            end
            viscE(iel) = max(viscE(iel), viscEq);
        end
    else
        R_max = max(abs(R));
        if abs(E_dev_max) < 1.0e-100
           viscE(iel) = 0.0;
        else
           viscE(iel) = cE*R_max / E_dev_max;
        end
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
