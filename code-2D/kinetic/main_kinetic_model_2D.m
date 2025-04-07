% Main file for 2D kinetic model
%   Psi_0 = Psi; Psi_1 = Phi \cdot \grad C
% No need to define domain: 
%   info of domain is in the file of initial data
%%
clc; clear;
mypara.project = 'kinetic'; mypara.init_ID = 1;
mypara.A_c = 20; mypara.eps = 1e-4;
time_duration = 20; time_step = 1e-2; time_initial = 0;
time_plotting = time_initial:0.1:time_initial+time_duration;
time_plotting_G = time_initial+time_duration;
Define_Function_2D; % myfunc
%% Initialization
% Load initial data
cd InitialData_2D/;
filename = strcat('data_2D_macro_ID_', num2str(mypara.init_ID), '.mat');
load(filename)
cd ..

% Neumann BC in v
domain.v_max = [5 5]; 
domain.v_min = -domain.v_max;
domain.v_number = [50 50];
domain.v_meshsize = (domain.v_max - domain.v_min) ./ domain.v_number;
domain.v1_vector = linspace(domain.v_min(1), domain.v_max(1), domain.v_number(1)+1);
domain.v2_vector = linspace(domain.v_min(2), domain.v_max(2), domain.v_number(2)+1);
% Initialize G
G1_initial = zeros([domain.x_number, domain.v_number+1]);
G2_initial = zeros([domain.x_number, domain.v_number+1]);

%% Time Evolution 
time_current = time_initial;
time_step_number = time_duration / time_step;
% Assign the initial value to be the current one
rho_current = rho_initial;
G1_current = G1_initial;
G2_current = G2_initial;
c_current = c_initial;
SAVE_DATA_2D;
% evolution
for kk_time = 1:time_step_number
    [rho_temp, c_temp, G1_temp, G2_temp] = OneStep_Kinetic2D(rho_current, c_current, G1_current, G2_current,...
        domain, time_step, myfunc, mypara);

    % update rho, c, G
    rho_current = rho_temp;
    c_current = c_temp;
    G1_current = G1_temp; 
    G2_current = G2_temp; 
    time_current = time_initial + time_step * kk_time;
    PLOT_DATA_2D;
    % save data and plot data   
    if min(abs(time_current-time_plotting)) < time_step/10 
        PLOT_DATA_2D;
        SAVE_DATA_2D;
    end
    
    if min(abs(time_current-time_plotting_G)) < time_step/10 
        SAVE_DATA_2D_G;
    end
end
