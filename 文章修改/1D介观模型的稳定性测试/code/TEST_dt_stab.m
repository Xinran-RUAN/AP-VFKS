% main file for testing the stability
% fix epsilon
mypara.eps = 1;
% test dt_list
% time_step_list = [ ...
%                   1e-4,2e-4,5e-4,...
%                   1e-3,2e-3,5e-3,...
%                   1e-2,2e-2,5e-2,...
%                   1e-1,2e-1,5e-1...
%                   ];
time_step_list = 1e-6;
N_dt_test = length(time_step_list);
for jj_dt=1:N_dt_test
    time_step = time_step_list(jj_dt);
    main_kinetic_model;
end