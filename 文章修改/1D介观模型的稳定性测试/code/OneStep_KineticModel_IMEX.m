% Periodic boundary: index 1 = index (N + 1)
%   where N = length(x)
% Index jj (starting from 1) ~ position in x
%       where x(kk) = x_min + domain.x_meshsize * kk, 
%             x_min = x(0) = x_1, (index 1)
%             x_max = x(N) = x_{N+1}. (index N+1)
%   rho, c:  x(0), x(1), ..., x(N-1); 
%   drho/domain.x_meshsize, dc/domain.x_meshsize, q(rho), q'(rho), Phi, g, tilde{g}:  x(1/2), x(3/2), ..., x(N-1/2); 
%   G1_flux, G2_flux, G_flux, Source_G: x(1/2), x(3/2), ..., x(N-1/2);
function[Rho_CurrentStep, C_CurrentStep, G_CurrentStep] = OneStep_KineticModel_IMEX(rho_current, c_current, G_current, domain, time_step, mypara, myfunc)
Nx = domain.x_number;
% compute phi and psi
Psi = myfunc.psi(domain.v_vector');
Phi = myfunc.phi(domain.v_vector');
% compute c, dc/domain.x_meshsize, rho, drho/domain.x_meshsize, q, dq/domain.x_meshsize, 
rho_current_extension = [rho_current, rho_current(1)];   % x(0), x(1), ..., x(N)
drho_dx = diff(rho_current_extension) / domain.x_meshsize;  
c_current_extension = [c_current, c_current(1)];         % x(0), x(1), ..., x(N)
dc_dx = diff(c_current_extension) / domain.x_meshsize;  
% compute q(rho) and q'(rho) at half grid points
rho_half = 0.5 * (rho_current_extension(1:end-1) + rho_current_extension(2:end)); 
q_half = myfunc.q_func(rho_half);    
dq_half = myfunc.dq_func(rho_half);  
% compute Phi = rho * q(rho) at half grid points
rho_left = rho_current;                          
rho_right = [rho_current(2:end), rho_current(1)]; 
PHI_half_explicit = myfunc.q_func(rho_right) .* rho_left .* (dc_dx >= 0) ...
            + myfunc.q_func(rho_left) .* rho_right .* (dc_dx < 0); 
% compute $G1_flux = v * \p_x(q(\rho)*g)$ and $G2_flux = v^2 * Psi * \p_x(rho*q'(rho)*\p_x(rho))$
qG = q_half .* G_current; 
qG_extension = [qG(:, end), qG, qG(:,1)];  
dqGdx_extension = (qG_extension(:, 2:end) - qG_extension(:, 1:end-1)) / domain.x_meshsize;   % x(0), x(1), ..., x(N)
G1_Flux = dqGdx_extension(:, 1:end-1) .* max(domain.v_vector', 0) ...
            + dqGdx_extension(:, 2:end) .* min(domain.v_vector', 0); % left or right derivative
        
rho_half_extension = [rho_half(end), rho_half, rho_half(1)]; % x(-1/2), x(1/2), x(3/2), ..., x(N-1/2), x(N+1/2)
Rho_dQ_dRhodx = rho_current_extension .* myfunc.dq_func(rho_current_extension) .* diff(rho_half_extension) / domain.x_meshsize; % x(0), x(1), ..., x(N)
G2_Flux = (domain.v_vector').^2 .* Psi .* diff(Rho_dQ_dRhodx) / domain.x_meshsize; 
% compute (I - Pi) (G1_flux + G2_flux)
G_Flux = G1_Flux + G2_Flux;
G_Flux_proj = integration_v_meshgrid(G_Flux, domain.v_meshsize) .* Psi; % projection -- sum over v
G_Flux_res = G_Flux - G_Flux_proj;   
% compute $r0*g*(1-rho/rho_m)_+$
Source_G = myfunc.f2_func(G_current, rho_half); 
%% determine dt
dt = time_step;
%% update tilde{G}
% compute tilde{S} without the term q*g/eps^2 
tS_remained = - domain.v_vector' .* Psi .* q_half .* drho_dx ...
    + Phi .* dc_dx .* PHI_half_explicit; 
% update g
tG = (dt * (-G_Flux_res / mypara.eps + tS_remained / mypara.eps ^ 2 + Source_G) + G_current) ./ (1 + dt / mypara.eps ^ 2 * q_half);    % x(1/2), x(3/2), ..., x(N-1/2)
%% update rho - implicit
Dh = integration_v_meshgrid((domain.v_vector.^2)' .* Psi, domain.v_meshsize);
% compute coefficients a
a = (dt * q_half) ./ (mypara.eps^2 + dt * q_half) .* Dh .* q_half;
a_coeff = a - Dh * rho_half .* dq_half;   % x(1/2), x(3/2), ..., x(N-1/2)
a_ext = [a_coeff(end), a_coeff];  % x(-1/2), x(1/2), x(3/2), ..., x(N-1/2)
% compute coefficients b
beta = integration_v_meshgrid(domain.v_vector' .* Phi, domain.v_meshsize);
b = (dt * q_half) ./ (mypara.eps^2 + dt * q_half) .* beta .* dc_dx; % x(1/2), x(3/2), ..., x(N-1/2)
b_ext = [b(end), b];  % x(-1/2), x(1/2), x(3/2), ..., x(N-1/2)
% construct matrix
Q_0 = [myfunc.q_func(rho_left(end)), myfunc.q_func(rho_left(1:end-1))];   % x(-1), x(0), x(1), ..., x(N-2)
Q_2 = myfunc.q_func(rho_right);  % x(1), ..., x(N-2), x(N-1), x(N)
m_0 = 1 / dt + (a_ext(1:end-1) + a_ext(2:end))' / domain.x_meshsize^2 ...
        + (max(b_ext(2:end), 0) .* Q_2 -  min(b_ext(1:end-1), 0) .* Q_0)' / domain.x_meshsize;
m_p1 = - a_ext(1:end-1)' / domain.x_meshsize^2 + (min(b_ext(1:end-1), 0) .* Q_0)' / domain.x_meshsize;
m_n1 = - a_coeff' / domain.x_meshsize^2 - (max(b_ext(2:end), 0) .* Q_2)' / domain.x_meshsize;
MAT_Rho = spdiags([m_n1, m_0, m_p1], [-1, 0, 1], Nx, Nx);
% Periodic B.C.
MAT_Rho(1, end) = m_n1(end);
MAT_Rho(end, 1) = m_p1(1);

% compute residuals "integral(v*diff(q*g), dv)" - index 1, 2, ..., N-1
q_tG = q_half .* tG; % x(1/2), x(3/2), ..., x(N-1/2)
q_tG_ext = [q_tG(:, end), q_tG];    
d_q_tG_dx = (q_tG_ext(:, 2:end) - q_tG_ext(:, 1:end-1)) / domain.x_meshsize;  
int_v_dqg = integration_v_meshgrid(domain.v_vector' .* d_q_tG_dx, domain.v_meshsize); % x(0), x(1), ..., x(N-1)
% compute residuals "diff(a*diff(rho))"
a_dRho = a .* drho_dx; % x(1/2), x(3/2), ..., x(N-1/2)
a_dRho_ext = [a_dRho(end), a_dRho]; %  x(-1/2), x(1/2), x(3/2), ..., x(N-1/2)
d_adrho = diff(a_dRho_ext) / domain.x_meshsize;
% compute source
r0_s = myfunc.f1_func(rho_current);
% compute chemotaxis term
b_PHI = b .* PHI_half_explicit; % 3/2, ..., N-1/2
b_PHI_ext = [b_PHI(end), b_PHI]; % 1/2, 3/2, ..., N-1/2
d_bPhi = diff(b_PHI_ext) / domain.x_meshsize;  % 1, 2, ..., N-1
% compute r
rhs = rho_current / dt - int_v_dqg - d_adrho + d_bPhi + r0_s;    % 1, ..., N-1
% update rho
Rho_CurrentStep_t = MAT_Rho \ rhs';
Rho_CurrentStep = Rho_CurrentStep_t';
%% update G
Rho_CurrentStep_extension = [Rho_CurrentStep, Rho_CurrentStep(1)];  % 1, 2, ..., N-1, N
drho_dx_CurrentStep = diff(Rho_CurrentStep_extension) / domain.x_meshsize;         % 3/2, ..., N-1/2

Rho_left_CurrentStep = Rho_CurrentStep;                                % 1, 2, ..., N-1
Rho_right_CurrentStep = [Rho_CurrentStep(2:end), Rho_CurrentStep(1)]; % 2, 3, ..., N-1, N
PHI_half_CurrentStep = myfunc.q_func(rho_right) .* Rho_left_CurrentStep .* (dc_dx >= 0) ...
                    + myfunc.q_func(rho_left) .* Rho_right_CurrentStep .* (dc_dx < 0);  % 3/2, ..., N-1/2


G_CurrentStep = tG + dt ./ (mypara.eps ^ 2 + dt * q_half) .* ( ...
    + Phi .* dc_dx .* (PHI_half_CurrentStep - PHI_half_explicit) ...
    - domain.v_vector' .* Psi .* q_half .* (drho_dx_CurrentStep - drho_dx));
%% update c
% zeta \p_t c = \Delta c + rho - c
% zeta is small since the chemoattractant diffuses faster...
zeta = 0;
% construct matrix
MAT_C = spdiags([-dt / domain.x_meshsize^2 * ones(Nx, 1), (dt * ( 1 + 2 / domain.x_meshsize^2) + zeta) * ones(Nx, 1), -dt / domain.x_meshsize^2 * ones(Nx, 1)], [-1, 0, 1], Nx, Nx);
% Periodic B.C.
MAT_C(1, end) = MAT_C(1, 2);
MAT_C(end, 1) = MAT_C(end, end-1);
% 
rhs = zeta * c_current' + dt * Rho_CurrentStep';
if issymmetric(MAT_C)
    opts.SYM = true;
    opts.POSDEF = true;
    C_CurrentStep_t = linsolve(full(MAT_C), rhs, opts);   
    % if MAT_C and rhs are both symmetric, we force C to be symmetric
    middle_ind = length(rhs) / 2 + 1;
    if rhs(end:-1:middle_ind+1) == rhs(2:middle_ind-1)
        C_CurrentStep_t(end:-1:middle_ind+1) = C_CurrentStep_t(2:middle_ind-1);
    end
else
    C_CurrentStep_t = MAT_C \ rhs;
end
C_CurrentStep = C_CurrentStep_t';
%%
%
% if issymmetric(MAT_Rho)
%     opts.SYM = true;
%     opts.POSDEF = true;
%     Rho_CurrentStep_t = linsolve(full(MAT_Rho), rhs', opts);
% 
%     % if MAT_Rho and rhs are both symmetric, we force Rho to be symmetric
%     middle_ind = length(rhs) / 2 + 1;
%     if rhs(end:-1:middle_ind+1) == rhs(2:middle_ind-1)
%         Rho_CurrentStep_t(end:-1:middle_ind+1) = Rho_CurrentStep_t(2:middle_ind-1);
%     end
% else
%     Rho_CurrentStep_t = MAT_Rho \ rhs';
% end