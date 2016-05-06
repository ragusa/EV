function l2norm = evaluate_entropy_jump_l2norm_ss(...
    u,mesh,phys,quadrature,ev,dof_handler)

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

% compute entropy jumps
n_face = n_dof;
jump = zeros(n_face,1);
integral = 0.0;
for f = 1:n_face % loop over interior faces
    if f == 1 || f == n_face
      % jumps on exterior faces considered zero
      continue
    else
      iL = f-1; % left  cell index
      iR = f;   % right cell index
    end
    uF = u(g(iL,2)); % solution on face
    dudx_L = dvdz' * u(g(iL,:)) / Jac(iL);
    dudx_R = dvdz' * u(g(iR,:)) / Jac(iR);
    dEdx_L = entropy_deriv(uF)*dudx_L(end);
    dEdx_R = entropy_deriv(uF)*dudx_R(1);
    jump(f) = abs(mu*(dEdx_L - dEdx_R));
end

% compute integral
integral = 0.0;
for iel = 1:nel
    % compute max jump
    faceL = g(iel,1);
    faceR = g(iel,2);
    jump_max = max(jump(faceL),jump(faceR));
    
    % add to integral
    integral = integral + sum(wq) * jump_max^2 * Jac(iel);
end

% compute square root to finish L^2 norm
l2norm = sqrt(integral);

