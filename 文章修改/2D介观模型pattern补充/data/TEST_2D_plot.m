% 读入2D数据并画图
%% 读入数据，记录程序是否可算出结果
% 若energy——record中出现NAN，表明程序崩溃
% eps_list = [1e-3,1e-2,1e-1,1];
eps_list = [2e-1,5e-1];
A = 20; % A = 20, 50
init_ID = 1; % init_ID = 1, 10
T_test = 20;
dt = 1e-3;

for jj = 1:length(eps_list)
    eps_plot = eps_list(jj);
    
    path_name = strcat('data_eps_',num2str(eps_plot, '%.0e'));
    filename_data = strcat('data_2D_kinetic_A_', num2str(A),...
                '_eps_', num2str(eps_plot, '%.0e'),...
                '_T_', num2str(T_test, '%.1f'), ...
                '_dt_1e-03_Nx_80_80_r_1e-01_init_', num2str(init_ID),'.mat');
    cd(path_name);
    load(filename_data);
    cd ..
    figure(jj);
    image(domain.x1_vector, domain.x2_vector, rho_current', 'CDataMapping','scaled'); 
    colorbar;
    caxis([0,1])
    title(strcat('$\chi_0$ = ', num2str(A), ', $\varepsilon$ =', num2str(eps_plot)),'Interpreter', 'latex');
    xlabel('$x_1$', 'Interpreter', 'latex');
    ylabel('$x_2$', 'Interpreter', 'latex');
    set(gca, 'FontSize', 20, 'LineWidth', 2);
    
    fig_name = strcat('Fig_2D_A_',num2str(A),'_eps_', num2str(eps_plot, '%.0e'),...
        '_T_', num2str(T_test, '%i'),...
        '_init_',num2str(init_ID),'.eps');
    saveas(gcf, fig_name, 'epsc');
end
