function uL = low_order_step(u_old,AL,ALtr_mod,ML,b,dt,theta,inc,periodic_BC)

rhs = (ML - (1-theta)*dt*AL)*u_old + dt*b;
if ~periodic_BC
    rhs(1) = inc;
end
uL = ALtr_mod \ rhs;

end