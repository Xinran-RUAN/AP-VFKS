% Macro model
% \p_t u = \p_x(D(u)\p_x u) - \p_x(A phi(u) \p_x c) + ru*(1-u/u_M)
% D(u) = q(u) - dq(u) * u
% phi(u) = q(u) * u
% The code is adapted to the Periodic BC 

%%
clc; clear;
mypara.A_c = 6; mypara.dataID = 1;
time_duration = 1; time_step = 1e-4; time_initial = 0;
time_plotting = time_initial:time_initial + time_duration;
Define_Function_macro;
%% initialization 
% SET_INITIAL_DATA;
cd InitialData/;
data_file = strcat('data_InitialTime_', num2str(time_initial, '%i'), ...
                '_Id_', num2str(mypara.dataID), '.mat');
load(data_file, 'domain', 'rho_initial', 'c_initial');
cd ..;

% No need to define domain since the domain is loaded from initial data
% Otherwise use "Define_Domain_macro;"

% Compute initial energy
[rho_ref, intH_ref] = Get_intH_ref(mypara, myfunc); 
%% Update
time_current = time_initial;
time_step_number = time_duration / time_step;
% Assign the initial value to be the current one
rho_current = rho_initial;
c_current = c_initial;
% initial energy
energy_record_time = time_initial:max(1e-2, time_step):time_initial + time_duration;
energy_current = Compute_Energy(rho_ref, intH_ref, rho_current, c_current, domain.x_meshsize);
energy_record = zeros(size(energy_record_time));
kk_energy = 1;
energy_record(kk_energy) = energy_current;
for kk_time = 1:time_step_number
    [rho_temp, c_temp] = OneStep_MacroModel_IMEX(rho_current, c_current,...
        domain, time_step, myfunc, mypara);
    rho_current = rho_temp;
    c_current = c_temp;
    time_current = time_current + time_step;
    if min(abs(time_current-energy_record_time)) < time_step / 2
        energy_current = Compute_Energy(rho_ref, intH_ref, rho_current, c_current, domain.x_meshsize);
        kk_energy = kk_energy + 1;
        energy_record(kk_energy) = energy_current;
        disp(time_current)
    end

    if min(abs(time_current - time_plotting)) < time_step/10
        PLOT_DATA;        
        SAVE_DATA_macro;
    end
end

