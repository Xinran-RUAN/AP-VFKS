% Plot density rho
figure(1)
image(domain.x1_vector, domain.x2_vector, rho_current', 'CDataMapping','scaled'); 
colorbar;
title(strcat('time =', num2str(time_current)));
xlabel('$x_1$', 'Interpreter', 'latex');
ylabel('$x_2$', 'Interpreter', 'latex');
set(gca, 'FontSize', 20, 'LineWidth', 2);
pause(0.1)