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
% Please check int_v psi(v) dv = 1 and D = int (v \cross v) psi(v) dv = Id(length(v), length(v))
% v1, v2 are row vectors
myfunc.psi = @(v1, v2) 1 / (2 * pi) * exp(-(v1' .* v1' + v2 .* v2 )/ 2); 
% Please check int_v v_j^2*phi(v) dv = A_c, int_v v_j*v_i*phi(v) dv = 0
myfunc.phi_1 = @(v1, v2) mypara.A_c / (2 * pi) * v1' .* exp(-(v1' .* v1' + v2 .* v2) / 2);
myfunc.phi_2 = @(v1, v2) mypara.A_c / (2 * pi) * v2 .* exp(-(v1' .* v1' + v2 .* v2) / 2);
% q(rho) and q'(rho)
myfunc.q_func = @(rho) max(1 - (rho / brho).^gm, 0) .* (rho > 0);
myfunc.dq_func = @(rho) -gm * (rho / brho).^(gm - 1) ./ brho .* (rho < brho) .* (rho > 0);
% proliferation function f
myfunc.f1_func = @(rho) r0 * rho .* max(1 - rho / mrho, 0).* (rho < brho) .* (rho > 0);
myfunc.f2_func = @(g, rho) r0 * g .* max(1 - rho / mrho, 0).* (rho < brho) .* (rho > 0);

%%
% wid = 1;
% myfunc.psi = @(v1, v2) 1 / (2 * pi) / wid^2 * exp(-(v1' .* v1' + v2 .* v2 )/ 2 / wid^2); 
% % Please check int_v v_j^2*phi(v) dv = A_c, int_v v_j*v_i*phi(v) dv = 0
% myfunc.phi_1 = @(v1, v2) mypara.A_c / (2 * pi) / wid^4 * v1' .* exp(-(v1' .* v1' + v2 .* v2) / 2 / wid^2);
% myfunc.phi_2 = @(v1, v2) mypara.A_c / (2 * pi) / wid^4 * v2 .* exp(-(v1' .* v1' + v2 .* v2) / 2 / wid^2);

