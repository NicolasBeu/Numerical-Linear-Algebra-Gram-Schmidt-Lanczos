% Problem 2: The Lanczos Iteration and Ritz Value Convergence Bounds
clear; clc; close all;

n = 100;
kmax = 20;

% Build Matrix A
A = zeros(n, n);
for i = 1:n
    for j = 1:n
        A(i,j) = min(i,j);
    end
end

true_eigs = sort(eig(A), 'descend');
lambda_1 = true_eigs(1);
lambda_n = true_eigs(end);

% Random Initial Vector
rng(0);
v = randn(n, 1);
v = v / norm(v);

% Run Lanczos Iteration for 20 iterations
[T, Q] = Lanczos_Desc(A, v, kmax);

% Compute Ritz values across iterations
ritz_values = zeros(kmax, kmax);
for j = 1:kmax
    Tj = T(1:j, 1:j);
    ritz_values(1:j, j) = sort(eig(Tj), 'descend');
end

% Error tracking for first Ritz value
mu1 = zeros(kmax, 1);
err1 = zeros(kmax, 1);
other_max = 455;
for j = 1:kmax
    mu1(j) = ritz_values(1, j);
    err1(j) = abs(mu1(j) - lambda_1);
end

% Upper Bound Calculation (Theorem 1)
rho = (lambda_1 - other_max) / (other_max - lambda_n);
ub = zeros(kmax, 1);
for j = 1:kmax
    term = (2 * rho^(j-2) / (1 + rho^(2*(j-2))))^2;
    ub(j) = (lambda_1 - lambda_n) * term;
end

% Error tracking for second Ritz value
mu2 = zeros(kmax, 1);
err2 = zeros(kmax, 1);
for j = 1:kmax
    if j > 1
        mu2(j) = ritz_values(2, j);
        err2(j) = abs(mu2(j) - lambda_1);
    else
        mu2(j) = NaN;
        err2(j) = NaN;
    end
end

% Plotting Convergence
figure('Name', 'Ritz Value Convergence', 'Position', [100, 100, 750, 500]);
semilogy(1:kmax, err1, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
semilogy(1:kmax, ub, 'r--', 'LineWidth', 1.5);
semilogy(1:kmax, err2, 's--', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Iteration (k)');
ylabel('Error in Ritz Value');
legend({'|\mu_1 - \lambda_1|', 'Upper Bound', '|\mu_2 - \lambda_1|'}, 'Location', 'northeast');
title('Ritz Values Convergence and Upper Bound Comparison');
grid on;

%% Supporting Functions for Problem 2

function [T, Q] = Lanczos_Desc(A, q, k)
    n = length(q);
    Q = zeros(n, k+1);
    T = zeros(k, k);
    
    Q(:,1) = q / norm(q);
    beta = 0;
    
    for j = 1:k
        z = A * Q(:,j);
        if j > 1
            z = z - beta * Q(:,j-1);
        end
        alpha = Q(:,j)' * z;
        z = z - alpha * Q(:,j);
        beta = norm(z);
        
        T(j,j) = alpha;
        if j < k
            T(j+1, j) = beta;
            T(j, j+1) = beta;
        end
        
        if beta == 0
            Q = Q(:, 1:j);
            T = T(1:j, 1:j);
            break;
        end
        Q(:, j+1) = z / beta;
    end
end
