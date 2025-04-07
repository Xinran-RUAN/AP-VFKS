  % \p_t u = \p_x(D(u^e)\p_x u^i) - \p_x(A phi(u) \p_x c) + ru*(1-u/u_M)
% The code is adapted to periodic BC
% index of unknowns: 1, ..., N-1 with 0 = N-1, 1 = N
% ex-implicit
function[rho_n, c_n] = OneStep_MacroModel_IMEX(rho_current, c_current, domain, time_step, myfunc, mypara)
Nx = length(rho_current);
beta = mypara.A_c;
%% update rho
% D_term - implicit
Rho_PreviousStep_extension = [rho_current(end), rho_current, rho_current(1)];    % 0, 1, ..., N-1, N
Q_PreviousStep_extension = myfunc.q_func(Rho_PreviousStep_extension);  % 0, 1, ..., N-1, N
Rho_half_PreviousStep_extension = 0.5 * (Rho_PreviousStep_extension(1:end-1) + Rho_PreviousStep_extension(2:end));  % 1/2, 3/2, ..., N-1/2
D_half_PreviousStep_extension = myfunc.d_func(Rho_half_PreviousStep_extension');   % 1/2, 3/2, ..., N-1/2

C_PreviousStep_extend = [c_current(end), c_current, c_current(1)];   % 0, 1, ..., N-1, N
b_dCdx_half_PreviousStep = beta * diff(C_PreviousStep_extend') / domain.x_meshsize;  % 1/2, 3/2, ..., N-1/2

a_n1 = -time_step / domain.x_meshsize^2 * [D_half_PreviousStep_extension(2:end-1); 0]...    % index of d: 3/2, ..., N-3/2, N-1/2(#)
    - time_step / domain.x_meshsize * max(b_dCdx_half_PreviousStep(2:end), 0) .* Q_PreviousStep_extension(3:end)';   % index of q: 2,3,..., N-1, N(#) % index of dx_C: 3/2, ..., N-3/2, N-1/2(#)
a_0 = ones(Nx, 1) + time_step / domain.x_meshsize^2 * (D_half_PreviousStep_extension(2:end) + D_half_PreviousStep_extension(1:end-1))...
    + time_step / domain.x_meshsize * (max(b_dCdx_half_PreviousStep(2:end), 0) .*  Q_PreviousStep_extension(3:end)' - min(b_dCdx_half_PreviousStep(1:end-1), 0) .*  Q_PreviousStep_extension(1:end-2)'); 
a_p1 = -time_step / domain.x_meshsize^2 * [0; D_half_PreviousStep_extension(2:end-1)]...    % index of d: 1/2(#), 3/2, ..., N-3/2
    + time_step / domain.x_meshsize * min(b_dCdx_half_PreviousStep(1:end-1), 0) .* Q_PreviousStep_extension(1:end-2)'; % index of q: 0(#), 1, ..., N-2 % index of dx_C: 1/2(#), 3/2, ..., N-3/2
MAT_D = spdiags([a_n1, a_0, a_p1], [-1, 0, 1], Nx, Nx);
MAT_D(1, end) = -time_step / domain.x_meshsize^2 * D_half_PreviousStep_extension(1) - time_step / domain.x_meshsize * max(b_dCdx_half_PreviousStep(1), 0) .* Q_PreviousStep_extension(2);  % 1/2
MAT_D(end, 1) = -time_step / domain.x_meshsize^2 * D_half_PreviousStep_extension(end) + time_step / domain.x_meshsize * min(b_dCdx_half_PreviousStep(end), 0) .* Q_PreviousStep_extension(end-1);    % N-1/2 (= 1/2)

% Source term
S_term = myfunc.f_func(rho_current');
% 
rho_n = (MAT_D \ (rho_current' + time_step * S_term))';
%% update c
zeta = mypara.zeta;
% construct matrix
MAT_C = spdiags([-time_step / domain.x_meshsize^2 * ones(Nx, 1), (time_step * ( 1 + 2 / domain.x_meshsize^2) + zeta) * ones(Nx, 1), -time_step / domain.x_meshsize^2 * ones(Nx, 1)], [-1, 0, 1], Nx, Nx);
% Periodic B.C.
MAT_C(1, end) = MAT_C(1, 2);
MAT_C(end, 1) = MAT_C(end, end-1);

c_n = MAT_C \ (zeta * c_current' + time_step * rho_n');
c_n = c_n';