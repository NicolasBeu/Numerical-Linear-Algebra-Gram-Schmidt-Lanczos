% Problem 1: Classical vs Modified Gram-Schmidt and Lanczos Orthonormalization
clear; clc; close all;

n = 100;
k = 6;
v = ones(n, 1);

% Build Matrix A where A(i,j) = min(i,j)
A = zeros(n);
for i = 1:n
    for j = 1:n
        A(i,j) = min(i,j);
    end
end

% Construct Krylov Subspace V
V = zeros(n, k);
V(:,1) = v;
for j = 2:k
    V(:,j) = A * V(:,j-1);
end

% Classical GS vs Modified GS execution
[Q_cgs, ~] = classicalGS(V);
[Q_mgs, ~] = modifiedGS(V);

orth_error_cgs = norm(Q_cgs' * Q_cgs - eye(k));
orth_error_mgs = norm(Q_mgs' * Q_mgs - eye(k));
cond_V = cond(V);

disp(['Condition number of V: ', num2str(cond_V)]);
disp(['Orthogonal Error Classical GS: ', num2str(orth_error_cgs)]);
disp(['Orthogonal Error Modified GS: ', num2str(orth_error_mgs)]);

% Lanczos Basis Comparison for dimension k = 6
[T_lanczos, Q_lanczos] = Lanczos_ex(A, v, k);
orth_error_lanczos = norm(Q_lanczos(:, 1:k)' * Q_lanczos(:, 1:k) - eye(k));
disp(['Orthogonal Error Lanczos: ', num2str(orth_error_lanczos)]);

%% Supporting Functions for Problem 1

function [Q, R] = classicalGS(V)
    [n, k] = size(V);
    Q = zeros(n, k);
    R = zeros(k, k);
    for j = 1:k
        q = V(:,j);
        for i = 1:j-1
            R(i,j) = Q(:,i)' * V(:,j);
            q = q - R(i,j) * Q(:,i);
        end
        R(j,j) = norm(q);
        if R(j,j) == 0
            error('linear dependent vectors in CGS');
        end
        Q(:,j) = q / R(j,j);
    end
end

function [Q, R] = modifiedGS(V)
    [n, k] = size(V);
    Q = V;
    R = zeros(k, k);
    for i = 1:k
        R(i,i) = norm(Q(:,i));
        if R(i,i) == 0
            error('linear dependent vectors in MGS');
        end
        Q(:,i) = Q(:,i) / R(i,i);
        for j = i+1:k
            R(i,j) = Q(:,i)' * Q(:,j);
            Q(:,j) = Q(:,j) - R(i,j) * Q(:,i);
        end
    end
end

function [T, Q] = Lanczos_ex(A, v, k)
    n = length(v);
    Q = zeros(n, k);
    qprev = v / norm(v);
    Q(:,1) = qprev;
    T = zeros(k, k);
    
    for ii = 1:k
        if ii == 1
            ui = A * qprev;
            T(1,1) = qprev' * ui;
            qi = ui - T(1,1) * qprev;
            if ii < k
                beta = norm(qi);
                T(ii+1, ii) = beta;
                T(ii, ii+1) = beta;
                qprev = qi / beta;
                Q(:, ii+1) = qprev;
            end
        else
            ui = A * Q(:,ii) - T(ii, ii-1) * Q(:,ii-1);
            T(ii, ii) = Q(:,ii)' * ui;
            qi = ui - T(ii, ii) * Q(:,ii);
            if ii < k
                beta = norm(qi);
                T(ii+1, ii) = beta;
                T(ii, ii+1) = beta;
                Q(:, ii+1) = qi / beta;
            end
        end
    end
end
