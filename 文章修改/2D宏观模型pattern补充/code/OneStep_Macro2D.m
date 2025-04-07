% IMEX scheme for
%   \p_t u = \nabla(D(u^{ex})\nabla u^{im}) - \nabla(A phi(u^{imex}) \nabla
%       c^{ex}) + f(u^{ex})
% The code is adapted to periodic BC
% index of unknowns: 1, ..., N-1 with 0 = N-1, 1 = N

function[rho_n, c_n] = OneStep_Macro2D(rho_current, c_current, domain, time_step, myfunc, mypara)
Nx = length(rho_current);
beta = mypara.A_c;
%% update rho
Nx_1 = domain.x_number(1);
Nx_2 = domain.x_number(2);
Nx_total = Nx_1 * Nx_2;
rho = rho_current;
c = c_current;

%% test D
Dh = 1; 
dt = time_step;
dx1 = domain.x_meshsize(1);
dx2 = domain.x_meshsize(2);

q_func = myfunc.q_func;
dq_func = myfunc.dq_func;

rho_1_h = 0.5 * (rho + rho([2:end,1], :)); % i: 1/2, ..., N-1/2
rho_2_h = 0.5 * (rho + rho(:, [2:end,1])); % j: 1/2, ..., N-1/2
q_1_h = q_func(rho_1_h); % 1/2, ..., N-1/2
q_2_h = q_func(rho_2_h); % 1/2, ..., N-1/2
dq_1_h = myfunc.dq_func(rho_1_h); % 1/2, ..., N-1/2
dq_2_h = myfunc.dq_func(rho_2_h); % 1/2, ..., N-1/2

D_r = Dh * q_1_h ... 
    - Dh * rho_1_h .* dq_1_h; 

D_u = Dh * q_2_h ...
    - Dh * rho_2_h .* dq_2_h; 

D_l = Dh * q_1_h([end, 1:end-1], :)...
    - Dh * rho_1_h([end, 1:end-1],:) .* dq_1_h([end, 1:end-1],:); 

D_d = Dh * q_2_h(:, [end, 1:end-1]) ...
    - Dh * rho_2_h(:,[end, 1:end-1]) .* dq_2_h(:,[end, 1:end-1]);

% rho_l = [rho(end, :); rho(1:end-1, :)]; % rho_{i-1,j}
% rho_r = [rho(2:end, :); rho(1,:)]; % rho_{i+1,j}
% rho_u = [rho(:, 2:end), rho(:, 1)]; % rho_{i,j+1}
% rho_d = [rho(:,end), rho(:, 1:end-1)]; % rho_{i,j-1}
% rho_lh = 0.5 * (rho + rho_l);
% rho_rh = 0.5 * (rho + rho_r);
% rho_uh = 0.5 * (rho + rho_u);
% rho_dh = 0.5 * (rho + rho_d);
% D_l = myfunc.q_func(rho_lh) - rho_lh .* myfunc.dq_func(rho_lh); 
% D_r = myfunc.q_func(rho_rh) - rho_rh .* myfunc.dq_func(rho_rh); 
% D_u = myfunc.q_func(rho_uh) - rho_uh .* myfunc.dq_func(rho_uh); 
% D_d = myfunc.q_func(rho_dh) - rho_dh .* myfunc.dq_func(rho_dh); 

%%
dx_c1_r = (c([2:end,1],:) - c) / dx1; % i: 1/2, 3/2, ..., N-1/2

dx_c2_u = (c(:,[2:end,1]) - c) / dx2; % j: 1/2, 3/2, ..., N-1/2

dx_c1_l = (c - c([end,1:end-1], :)) / dx1; % i: -1/2, 1/2, ..., N-3/2

dx_c2_d = (c - c(:, [end,1:end-1])) / dx2; % j: -1/2, 1/2, ..., N-3/2

K_r =  beta * dt / dx1 * min(dx_c1_r, 0) .* q_func(rho);

K_u =  beta * dt / dx2 * min(dx_c2_u, 0) .* q_func(rho);

K_l = -beta * dt / dx1 * max(dx_c1_l, 0) .* q_func(rho);

K_d = -beta * dt / dx2 * max(dx_c2_d, 0) .* q_func(rho);

K_0 = beta * dt / dx1 * max(dx_c1_r, 0) .* q_func(rho([2:end,1], :)) ...
    + beta * dt / dx2 * max(dx_c2_u, 0) .* q_func(rho(:, [2:end,1]))...
    - beta * dt / dx1 * min(dx_c1_l, 0) .* q_func(rho([end,1:end-1], :)) ...
    - beta * dt / dx2 * min(dx_c2_d, 0) .* q_func(rho(:,[end,1:end-1]));


% c1_ex = [c(end, :); c; c(1, :)];
% dx_c1 = (c1_ex(2:end, :) - c1_ex(1:end-1, :)) / domain.x_meshsize(1);
% K_r = beta * time_step / domain.x_meshsize(1) * min(dx_c1(2:end, :), 0) .* myfunc.q_func(rho);
% K_l = -beta * time_step / domain.x_meshsize(1) * max(dx_c1(1:end-1, :), 0) .* myfunc.q_func(rho);
% c2_ex = [c(:, end), c, c(:, 1)];
% dx_c2 = (c2_ex(:, 2:end) - c2_ex(:, 1:end-1)) / domain.x_meshsize(2);
% K_u = beta * time_step / domain.x_meshsize(2) * min(dx_c2(:, 2:end), 0) .* myfunc.q_func(rho);
% K_d = -beta * time_step / domain.x_meshsize(2) * max(dx_c2(:, 1:end-1), 0) .* myfunc.q_func(rho);
% 
% K_0 = beta * time_step / domain.x_meshsize(1) * max(dx_c1(2:end, :), 0) .* myfunc.q_func(rho_r) ...
%     -beta * time_step / domain.x_meshsize(1) * min(dx_c1(1:end-1, :), 0) .* myfunc.q_func(rho_l) ...
%     +beta * time_step / domain.x_meshsize(2) * max(dx_c2(:, 2:end), 0) .* myfunc.q_func(rho_u)...
%     -beta * time_step / domain.x_meshsize(2) * min(dx_c2(:, 1:end-1), 0) .* myfunc.q_func(rho_d);

%%
% construct matrix of coefficients
A_l = -time_step / domain.x_meshsize(1)^2 * D_l + K_l; 
A_r = -time_step / domain.x_meshsize(1)^2 * D_r + K_r; 
A_d = -time_step / domain.x_meshsize(2)^2 * D_d + K_d;
A_u = -time_step / domain.x_meshsize(2)^2 * D_u + K_u; 
A_0 = 1 + time_step / domain.x_meshsize(1)^2 * (D_l + D_r) ...
        + time_step / domain.x_meshsize(2)^2 * (D_d + D_u) ...
        + K_0;

% construct MAT_D   
a_l = my_reshape(A_l, [Nx_1, Nx_2]);
a_ls = [a_l(2:end); 0]; % shifted for sparse
a_r = my_reshape(A_r, [Nx_1, Nx_2]);
a_rs = [0; a_r(1:end-1)];
a_d = my_reshape(A_d, [Nx_1, Nx_2]);
a_ds = [a_d(Nx_1+1:end); zeros(Nx_1,1)];
a_u = my_reshape(A_u, [Nx_1, Nx_2]);
a_us = [zeros(Nx_1,1); a_u(1:end-Nx_1)];
a_0 = my_reshape(A_0, [Nx_1, Nx_2]);
 
MAT_D = spdiags([a_ds, a_ls, a_0, a_rs, a_us], [-Nx_1, -1, 0, 1, Nx_1], Nx_total, Nx_total);
for jj = 2:Nx_1
    MAT_D(1+(jj-1)*Nx_1, (jj-1)*Nx_1) = 0;
end
for jj = 1:Nx_1-1
    MAT_D(jj*Nx_1, jj*Nx_1+1) = 0;
end
% periodic BC
for jj = 1:Nx_1
    MAT_D(1+(jj-1)*Nx_1, jj*Nx_1) = a_l(1+(jj-1)*Nx_1);
    MAT_D(jj*Nx_1, 1+(jj-1)*Nx_1) = a_r(jj*Nx_1);
    MAT_D(jj, end-Nx_1+jj) = a_d(jj);
    MAT_D(end-Nx_1+jj, jj) = a_u(end-Nx_1+jj);
end

% Source term
S_term = myfunc.f1_func(rho_current);
S_vec = reshape(S_term, [Nx_total, 1]);
% solve rho
rho_vec = my_reshape(rho_current, [Nx_1, Nx_2]);
rho_vec_n = MAT_D \ (rho_vec + time_step * S_vec);
rho_n = reshape(rho_vec_n, [Nx_1, Nx_2]);
%% update c
zeta = mypara.zeta;
% construct matrix
C_l = -time_step / domain.x_meshsize(1)^2 * ones(size(c));
C_r = -time_step / domain.x_meshsize(1)^2 * ones(size(c));
C_d = -time_step / domain.x_meshsize(2)^2 * ones(size(c));
C_u = -time_step / domain.x_meshsize(2)^2 * ones(size(c));
C_0 = zeta + time_step - (C_l + C_r + C_u + C_d);

c_l = my_reshape(C_l, [Nx_1, Nx_2]);
c_ls = [c_l(2:end); 0];
c_r = my_reshape(C_r, [Nx_1, Nx_2]);
c_rs = [0; c_r(1:end-1)];
c_d = my_reshape(C_d, [Nx_1, Nx_2]);
c_ds = [c_d(Nx_1+1:end); zeros(Nx_1,1)];
c_u = my_reshape(C_u, [Nx_1, Nx_2]);
c_us = [zeros(Nx_1,1); c_u(1:end-Nx_1)];
c_0 = my_reshape(C_0, [Nx_1, Nx_2]);

MAT_C = spdiags([c_ds, c_ls, c_0, c_rs, c_us], [-Nx_1, -1, 0, 1, Nx_1], Nx_total, Nx_total);
for jj = 2:Nx_1
    MAT_C(1+(jj-1)*Nx_1, (jj-1)*Nx_1) = 0;
end
for jj = 1:Nx_1-1
    MAT_C(jj*Nx_1, jj*Nx_1+1) = 0;
end
% Periodic B.C.
for jj = 1:Nx_1
    MAT_C(1+(jj-1)*Nx_1, jj*Nx_1) = c_l(1+(jj-1)*Nx_1);
    MAT_C(jj*Nx_1, 1+(jj-1)*Nx_1) = c_r(jj*Nx_1);
    MAT_C(jj, end-Nx_1+jj) = c_d(jj);
    MAT_C(end-Nx_1+jj, jj) = c_u(end-Nx_1+jj);
end
% solve c
c_vec = reshape(c, [Nx_total, 1]);
c_vec_n = MAT_C \ (zeta * c_vec + time_step * rho_vec_n);
c_n = reshape(c_vec_n, [Nx_1, Nx_2]);
