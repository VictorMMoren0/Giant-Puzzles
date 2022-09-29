clear all
close all
clc

%% Comments
%
% This code initially generates a random matrix of size n x n with entries
% 1 or -1, then at each iteration it generates 'm' random vectors of length
% n with entries 1 or -1 and it is applied the product of the matrix with
% the vector. After that, the corresponding vector to which the final
% vector (obtained by the product) has more zeros (more 'orthogonality'
% with the lines of the matrix) is substituted with the line that has the 
% greater inner product (in absolute value, less 'orthogonal'). If the
% absolute value of the determinant of the new matrix is greater than the 
% previous one, this is our new matrix and the code proceed to the next
% iteration, otherwise the matrix keeps the same.
% 
% Some additional parameters are given:
% - n_vec: if wanted, it is possible to create 'n_vec' vectors at each
% iteration to substitute 'n_vec' lines in the matrix.
% 
% - c_max: Since this process does not garantee the convergence to the
% maximum determinant, if the solution get stuck at some value (which can
% or cannot be the exact solution), it saves the matrix 'final_A', and
% starts with a new matrix.

%% General Parameters
n = 15; % Matrix dimensions
m = 10; % Number of random vectors at each iteration
max_itera = 1e2; % Maximum number of iterations

%% Additional Parameters
n_vec = 1;
c_max = 1e6;

%% Main
final_A = zeros(n);
itera = 0;

for itera = 1:max_itera
    A = 2 * randi([0,1],n) - 1; % Initial Matrix
    c = 0;
    
    while c < c_max
        if abs(det(A)) > abs(det(final_A))
            final_A = A;
        end

        max_zeros = 0;
        for i = 1:m
            random_vec = 2 * randi([0,1],[n_vec,n]) - 1;
            inner_prod_vec = A * random_vec';
            if sum(inner_prod_vec == 0) > max_zeros
                max_zeros = sum(inner_prod_vec == 0);
                random_vec_min = random_vec;
            end
        end

        [~,k] = max(abs(inner_prod_vec));
        new_A = A;
        new_A(k,:) = random_vec;

        c = c + 1;
        
        if abs(det(new_A)) > abs(det(A))
            A = new_A;
            [det(A) det(final_A)] % Show the current max determinant and the global max
            c = 0; % Restart the counting
        end
    end
end

if sign(det(final_A)) == -1
    % If the determinant is negative, change the first 2 lines so that
    % it becomes positive
    final_A([1,2],:) = final_A([2,1],:)
    factor(int(det(final_A)))
end
