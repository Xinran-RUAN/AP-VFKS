% Main file for 2D macro model
%%
clc; clear;
mypara.project = 'macro'; 
mypara.init_ID = 1; % ID = 1: 0.5 10: 0.1
mypara.A_c = 50; % Ac = 20, 50
time_duration = 20; time_step = 1e-3; time_initial = 0;
% time_plotting = time_initial:0.1:time_initial+time_duration;
time_plotting = time_initial:time_initial+time_duration;
Define_Function_2D; % myfunc
%% initialization 
% % Set initial data
% SET_INITIAL_DATA_2D; % domain, rho_init and c_init

% load initial data
cd InitialData_2D/
filename = strcat('data_2D_macro_ID_', num2str(mypara.init_ID), '.mat');
load(filename)
cd ..

% % for test
% x1v = domain.x1_vector;
% x2v = domain.x2_vector;
% rho_initial = exp(-0.1 * (x1v' .* x1v' + x2v .* x2v));
% c_initial = exp(-(x1v' .* x1v' + x2v .* x2v));
%% Update
time_current = time_initial;
time_step_number = time_duration / time_step;
% Assign the initial value to be the current one
rho_current = rho_initial;
c_current = c_initial;

SAVE_DATA_2D; % save initial data
%  evolution
for kk_time = 1:time_step_number
    [rho_temp, c_temp] = OneStep_Macro2D(rho_current, c_current,...
        domain, time_step, myfunc, mypara);
    rho_current = rho_temp;
    c_current = c_temp;
    time_current = time_initial + time_step * kk_time;

    if min(abs(time_current - time_plotting)) < time_step/10
        PLOT_DATA_2D;        
        SAVE_DATA_2D;
    end
end

%% change plot settings
colorbar;
caxis([0,1])
title(strcat('$\chi_0$ = ', num2str(mypara.A_c), ', macro'), 'Interpreter', 'latex');
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
set(gca, 'FontSize', 20, 'LineWidth', 2);
