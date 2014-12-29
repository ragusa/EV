function [uH,DH] = high_order_step(u_older,u_old,viscL,dx,x,omega,sigma,...
    source,inc,dt,v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,...
    A,b,MC,theta,time_step,high_order_scheme,periodic_BC)

% size of system
n_dof = length(u_old);

% compute high-order diffusion matrix
if (high_order_scheme == 2) % Entropy viscosity
    if (time_step == 1)
        DH = compute_high_order_diffusion_matrix(viscL,viscL,dx,periodic_BC);
    else
        viscE = compute_entropy_viscosity(...
            u_older,u_old,x,omega,sigma,source,dt,...
            v,dvdz,zq,wq,Jac,cE,cJ,entropy,entropy_deriv,periodic_BC);
        DH = compute_high_order_diffusion_matrix(viscE,viscL,dx,periodic_BC);
    end
else
    DH = spalloc(n_dof,n_dof,3*n_dof);
end

% compute high-order matrices
AH = A + DH;
AHtr = MC + dt*theta*AH;
AHtr_mod = AHtr;
if ~periodic_BC
    AHtr_mod(1,:)=0; AHtr_mod(1,1)=1;
end

% compute high-order solution
rhs = (MC - (1-theta)*dt*AH)*u_old + dt*b;
if ~periodic_BC
    rhs(1) = inc;
end
uH = AHtr_mod \ rhs;

end