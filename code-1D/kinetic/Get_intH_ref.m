% To numerically compute int_H(rho), the antiderivative of H(rho)
%   where H(rho) = D0/beta*ln(rho/q(rho))
% We set int_H(0.5) = 0, D0 = 0.5, beta = A_c
% Make sure rho is between 0 and 1
% jj: index in the vector, starting from 1
%   H:      range of rho: (1/2:1:N-1/2)*drho;  index jj ~ H((jj-1/2)*drho)
%   int_H:  range of rho: (0:1:N)*drho;        index jj ~ int_H((jj-1)*drho)
% We finally get a vector intH_ref with values int_H(rho = (0:1:N)*drho)
%   To get int_H(rho) numerically for a general rho, 
%       we combine intH_ref and linear interpolatopn.
%%
function[rho_ref, intH_ref] = Get_intH_ref(mypara, myfunc)
drho_ref = 1e-3;
rho_ref = 0:drho_ref:1;
rho_h_ref = (rho_ref(1:end-1) + rho_ref(2:end)) / 2;
qrho_h_ref = myfunc.q_func(rho_h_ref);
% Compute H(rho)
D0 = 0.5; bt = mypara.A_c; 
H_ref_half = D0 / bt * log(rho_h_ref ./ qrho_h_ref);    
% Set int_H(0.5)=0
intH_ref = zeros(size(rho_ref));                        
N_ref = (length(rho_ref) - 1) / 2;  % 2 * N_ref * drho = 1
mid = N_ref + 1;                    % (mid-1) * drho = 0.5
intH_ref(mid) = 0;
% Looping over (0,1) to compute int_H
for kk = 1:N_ref
    % int_H((mid+kk-1)*drho) = int_H((mid+kk-2)*drho) + H((mid+kk-3/2)*drho) * drho
    intH_ref(mid + kk) = intH_ref(mid + kk - 1) + H_ref_half(mid + kk - 1) * drho_ref;
    % int_H((mid-kk-1)*drho) = int_H((mid-kk)*drho) + H((mid-kk-1/2)*drho) * drho
    intH_ref(mid - kk) = intH_ref(mid - kk + 1) - H_ref_half(mid - kk) * drho_ref;
end
