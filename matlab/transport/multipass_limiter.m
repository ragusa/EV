function Flim = multipass_limiter(F,Qplus,Qminus,opts,fct_opts)

% size
n = size(F,1);

% unpack options
tol = fct_opts.multipass_tol;
dirichlet_limiting_coefficient = fct_opts.dirichlet_limiting_coefficient;

% compute total possible antidiffusion
total_antidiffusion = 0.5*sum(sum(abs(F),1));

% initialize
Flim = zeros(n,1);
F_remainder = F;
iter = 0;

% loop until very little antidiffusion accepted
while (true)
    % increment iteration number
    iter = iter + 1;

    % compute limiting coefficients and antidiffusion source
    [Flim_new,L] = fct_opts.limiter(F_remainder,Qplus,Qminus,Flim,...
        opts,dirichlet_limiting_coefficient);
    
    % update total antidiffusion source
    Flim = Flim + Flim_new;
    
    % compute accepted antidiffusion in this iteration
    F_accepted = F_remainder.*L;
    accepted_antidiffusion = 0.5*sum(sum(abs(F_accepted),1));
    percent_accepted = accepted_antidiffusion / total_antidiffusion * 100;

    % report
    fprintf('          limiter pass %i: %f%% antidiffusion accepted\n',...
      iter,percent_accepted);
    
    % check accepted antidiffusion against tolerance
    if (percent_accepted < tol)
        break;
    end
    
    % update remainder antidiffusive fluxes
    F_remainder = F_remainder - F_accepted;
end

end
