function uFCT = compute_uFCT(uL,L,flux)

% get low-order steady-state matrix and modify for Dirichlet BC
AL = importdata('AL.txt');
AL(1,:) = 0; AL(:,1) = 0; AL(1,1) = 1.0;

% get number of edges
n_dofs = length(uL);
n_edges = n_dofs - 1;

% compute limited fluxes
limited_fluxes = zeros(n_dofs,1);
for iedge = 1:n_edges
    limited_fluxes(iedge) = limited_fluxes(iedge) + L(iedge)*flux(iedge);
    limited_fluxes(iedge+1) = limited_fluxes(iedge+1) - L(iedge)*flux(iedge);
end
% modify for Dirichlet BC
limited_fluxes(1) = 0.0;

% compute solution correction
du = AL\limited_fluxes;

% compute FCT solution
uFCT = uL + du;

end