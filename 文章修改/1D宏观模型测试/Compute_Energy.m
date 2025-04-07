% compute energy
function[energy_current] = Compute_Energy(rho_ref, intH_ref, rho_current, c_current, x_meshsize)
intH_current = interp1(rho_ref, intH_ref, rho_current);
energy_current = sum(intH_current) * x_meshsize - 0.5 * sum(rho_current .* c_current) * x_meshsize;

