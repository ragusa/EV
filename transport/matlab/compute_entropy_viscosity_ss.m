function viscE = compute_entropy_viscosity_ss(...
    u,x,mu,sigma,source,v,dvdz,wq,Jac,cE,cJ,entropy,entropy_deriv)

n = length(u); % number of dofs
nel = n-1;         % number of elements
g = [linspace(1,nel,nel)' linspace(2,nel+1,nel)']; % connectivity array

% compute domain average of entropy
E_integral = 0;
for iel = 1:nel
    u_local = v' * u(g(iel,:)); % solution evaluated at quad. points
    Eq = entropy(u_local);
    E_integral = E_integral + dot(wq,Eq)*Jac(iel);
end
L = x(end)-x(1);        % length of domain
E_avg = E_integral / L; % domain average of entropy

% compute max entropy deviation in domain for normalization constant
E_dev_max = 0;
for iel = 1:nel
    u_local = v' * u(g(iel,:)); % solution evaluated at quad. points
    Eq = entropy(u_local);
    E_dev_max = max(E_dev_max, max(abs(Eq-E_avg)));
end

% compute entropy jumps
n_faces=n;
jump = zeros(n_faces,1);

for f = 1:n_faces % loop over interior faces
    if f==1 || f==n_faces
        continue
    else % interior face when Dirchlet
        iL = f-1;   % left  cell index
        iR = f; % right cell index
    end
    uF = u(g(iL,2)); % solution on face
    dudx_L = dvdz' * u(g(iL,:)) / Jac(iL);
    dudx_R = dvdz' * u(g(iR,:)) / Jac(iR);
    dEdx_L = entropy_deriv(uF)*dudx_L(end);
    dEdx_R = entropy_deriv(uF)*dudx_R(1);
    jump(f) = abs(mu*(dEdx_L - dEdx_R));
end

% compute entropy viscosity
viscE = zeros(nel,1);
for iel = 1:nel
    % evaluate new and old solution at quadrature points
    u_local     = v'    * u(g(iel,:));
    dudx_local  = dvdz' * u(g(iel,:)) / Jac(iel);
    % evaluate entropy and derivative
    E    = entropy(u_local);
    dEdu = entropy_deriv(u_local);
    dEdx = dEdu.*dudx_local;
    
    % compute maximum entropy residual at quadrature points
    R = mu*dEdx + dEdu.*(sigma(iel)*u_local - source(iel));
    R_max = max(abs(R));
    
    % compute max jump
    iel_next = iel+1;
    jump_max = max(jump(iel),jump(iel_next));
    
    % compute entropy viscosity
    viscE(iel) = (cE*R_max + cJ*jump_max)/ E_dev_max;
end