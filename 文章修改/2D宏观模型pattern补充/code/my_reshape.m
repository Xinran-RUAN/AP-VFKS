% To convert matrix to a vector
%   order along each direction
% In 2D: 
%   A_ij = a(x_i, x_j)
%   A_vec(i + n_i * (j-1)) = A_ij
function[A_vec] = my_reshape(A, A_size)
new_size = [prod(A_size), 1];
A_vec = reshape(A, new_size);