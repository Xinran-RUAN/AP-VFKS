% Main file for kinetic model
%   Psi_0 = Psi; Psi_1 = Phi * dx_C
% Parameters to be provided:
%   mypara.eps; mypara.A_c; N_data_initial (run main_code_macro.m first)
%   N_data_initial: Id for a random initial data
% No need to define domain: 
%   info of domain is in the file of initial data
%% time_step and eps not given here

%%
mypara.A_c = 20;  
mypara.dataID = 1;
time_duration = 23; 
time_initial = 17;
time_plotting = time_initial:time_initial + time_duration;
Define_Function;
%% Initialization
% SET_INITIAL_DATA;
% Load initial data from datafiles
cd InitialData/;
% data_file = strcat('data_InitialTime_', num2str(time_initial, '%i'), ...
%                 '_Id_', num2str(mypara.dataID), '.mat');
% load(data_file, 'domain', 'rho_initial', 'G_initial', 'c_initial');
load('data_kinetic_A_20_eps_1e+00_dt_1e-03_T_17_Nx_400_r_1e-01_init_1.mat');
rho_initial = rho_current;
G_initial = G_current;
c_initial = c_current;
cd ..;

% Compute initial energy
[rho_ref, intH_ref] = Get_intH_ref(mypara, myfunc); 
%% Time Evolution 
time_current = time_initial;
time_step_number = time_duration / time_step;
% Assign the initial value to be the current one
rho_current = rho_initial;
G_current = G_initial;
c_current = c_initial;
% initial energy
energy_record_time = time_initial:max(1e-2, time_step):time_initial + time_duration;
energy_current = Compute_Energy(rho_ref, intH_ref, rho_current, c_current, domain.x_meshsize);
energy_record = zeros(size(energy_record_time));
kk_energy = 1;
energy_record(kk_energy) = energy_current;

for kk_time = 1:time_step_number
    [rho_temp, c_temp, G_temp] = OneStep_KineticModel_IMEX(rho_current, c_current, G_current,...
        domain, time_step, mypara, myfunc);

    % update rho, c, G
    rho_current = rho_temp;
    c_current = c_temp;
    G_current = G_temp;   
    time_current = time_current + time_step;
    
    % compute energy
    if min(abs(time_current-energy_record_time)) < time_step / 2
        energy_current = Compute_Energy(rho_ref, intH_ref, rho_current, c_current, domain.x_meshsize);
        kk_energy = kk_energy + 1;
        energy_record(kk_energy) = energy_current;
        disp(time_current)
    end
    
    % save data and plot data   
    if min(abs(time_current-time_plotting)) < time_step / 2  
        PLOT_DATA;
        SAVE_DATA;
    end
end
