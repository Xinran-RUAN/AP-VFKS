% Periodic BC in x
domain.x_max = [10 10]; 
domain.x_min = -domain.x_max;
domain.x_number = [80, 80];
domain.x_meshsize = (domain.x_max - domain.x_min) ./ domain.x_number;
domain.x1_vector = linspace(domain.x_min(1), domain.x_max(1) - domain.x_meshsize(1), domain.x_number(1));
domain.x2_vector = linspace(domain.x_min(2), domain.x_max(2) - domain.x_meshsize(2), domain.x_number(2));


