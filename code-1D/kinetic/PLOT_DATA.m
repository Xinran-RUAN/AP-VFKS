% Plot rho and c
figure(1)
plot(domain.x_vector, rho_current, 'LineWidth', 2, 'Color', 'b');
hold on;
plot(domain.x_vector, c_current, '--', 'LineWidth', 2);
hold off   
title(strcat('time =', num2str(time_current)));
xlabel('$x$', 'Interpreter', 'latex');
set(gca, 'FontSize', 20, 'LineWidth', 2);
legend('$\rho(x)$', '$c(x)$', 'Interpreter', 'latex', 'Location', 'best');
pause(0.1)