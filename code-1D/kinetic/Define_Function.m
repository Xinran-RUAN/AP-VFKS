%% general parameters
mypara.gamma = 1; 
mypara.r0 = 0.1; 
mypara.zeta = 0;    
% Please check that bar_rho >= rho_max
mypara.bar_rho = 1;
mypara.rho_max = 0.5;
% localize parameters
gm = mypara.gamma; brho = mypara.bar_rho;
r0 = mypara.r0; mrho = mypara.rho_max;
%% functions
% Please check int_v psi(v) dv = 1 and D = 1
myfunc.psi = @(v) 1 / sqrt(2 * pi) * exp(-v .* v / 2); 
% Please check int_v phi(v) dv = A_c
myfunc.phi = @(v) 1 / sqrt(2 * pi) * v .* exp(-v .* v / 2) * mypara.A_c;
% q(rho) and q'(rho)
myfunc.q_func = @(rho) max(1 - (rho / brho).^gm, 0) .* (rho > 0);
myfunc.dq_func = @(rho) -gm * (rho / brho).^(gm - 1) ./ brho .* (rho < brho) .* (rho > 0);
% proliferation function f
myfunc.f1_func = @(rho) r0 * rho .* max(1 - rho / mrho, 0).* (rho < brho) .* (rho > 0);
myfunc.f2_func = @(g, rho) r0 * g .* max(1 - rho / mrho, 0).* (rho < brho) .* (rho > 0);