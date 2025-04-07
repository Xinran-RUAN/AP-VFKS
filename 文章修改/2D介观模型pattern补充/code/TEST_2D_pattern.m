% main file for testing patterns
% fix parameters
mypara.eps = 1;
time_step = 1e-3;
%%
InitID_list = [1,10];
Ac_list = [20, 50];
N_ID_test = length(InitID_list);
N_Ac_test = length(Ac_list);
for jj_ID=1:N_ID_test
    for jj_Ac = 1:N_Ac_test
        mypara.init_ID = InitID_list(jj_ID);
        mypara.A_c = Ac_list(jj_Ac);
        main_kinetic_model_2D;
    end
end