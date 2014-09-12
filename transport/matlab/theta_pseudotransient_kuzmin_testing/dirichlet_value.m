function val = dirichlet_value( n,dt,inc,theta )

time_n  =dt*(n-1);
time_np1=dt*(n  );
ramp=0.1;
val=   theta *min(inc,ramp*time_n)+...
    (1-theta)*min(inc,ramp*time_np1);

end

