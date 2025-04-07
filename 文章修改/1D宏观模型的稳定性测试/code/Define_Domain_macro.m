% Periodic BC in x
domain.x_max = 20; 
domain.x_number = 400; 
domain.x_min = -domain.x_max;
domain.x_meshsize = (domain.x_max - domain.x_min) / domain.x_number;
% domain.x_vector = linspace(domain.x_min, domain.x_max - domain.x_meshsize, domain.x_number);
mid = domain.x_number / 2 +1;
x = zeros(1, domain.x_number);
x(mid) = 0;
for kk = 1:mid-2
    x(mid+kk) = x(mid) + kk * domain.x_meshsize;
    x(mid-kk) = x(mid) - kk * domain.x_meshsize;
end
x(1) = domain.x_min;    % symmetry
domain.x_vector = x;
