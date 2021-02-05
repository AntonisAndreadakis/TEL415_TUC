%TEL415 - ERGASIA 3
%ARI8MOS OMADAS: LAB41544983
%Matsatsos Ioannis 2013030148
%Andreadakis Antonis 2013030059
clc; clear all; close all;

%% ASKHSH 2
load('real_data.mat')
fprintf('\n');
fprintf('\tWe will use a matrix decomposition method for real data.\n');
fprintf('\tWhich method would you like to use?\n');
fprintf('\tPress 1 for eigendecomposition or 2 for Cholesky:\n');
method = input(' ');
N = 10000; % N is number of independent vectors
X = real_colored_noise(m, C, N, method);
real_mean = mean(X,2);
fprintf('\tmean vector:\n');
disp(real_mean);
real_C = cov(X');
fprintf('\tcovariance matrix:\n');
disp(real_C);

load('complex_data.mat')
fprintf('\n');
fprintf('\tWe will use a matrix decomposition method for complex data.\n');
fprintf('\tWhich method would you like to use?\n');
fprintf('\tPress 1 for eigendecomposition or 2 for Cholesky:\n');
method = input(' ');
N = 10000; % N is number of independent vectors
X2 = complex_colored_noise(m,C,N,method);
complex_mean = mean(X2,2);
fprintf('\tmean vector:\n');
disp(complex_mean);
complex_C = cov(X2');
fprintf('\tcovariance matrix:\n');
disp(complex_C);


function X = real_colored_noise(m, C, N, method)
% m is vector n-by-1 and defines var of Gaussian vector
% C is n-by-n covariance matrix
n = length(C);
 % N is number of independent vectors
Y = randn(n,N); %Gaussian vector n-by-1
% choose method for decomposition
if(method == 1) % eigen decomposition
    fprintf('\tYou chose eigendecomposition.\n');
    [Q, L] = eig(C);
    F = Q*((L)^(1/2));
    res = F;
elseif(method == 2) % Cholesky decomposition
        fprintf('\tYou chose Cholesky decomposition.\n');
        L = chol(C);
        res = L;
        else
        fprintf('Invalid method. Choose 1 for eig or 2 for chol.\n');
 end
        Z = res*Y;
        X = Z + m.*ones(n,N);
        fprintf('\tThe size of matrix X is %dx%d.\n', n, N);
end

function X = complex_colored_noise(m,C,N,method)
% m is vector n-by-1 and defines var of Gaussian vector
% C is n-by-n covariance matrix
n = length(C);
 % N is number of independent vectors
Y = randn(n,N); %Gaussian vector n-by-1
% choose method for decomposition
if(method == 1) % eigen decomposition
    fprintf('\tYou chose eigendecomposition.\n');
    [Q, L] = eig(C);
    F = Q*((L)^(1/2));
    res = F;
elseif(method == 2) % Cholesky decomposition
        fprintf('\tYou chose Cholesky decomposition.\n');
        L = chol(C);
        res = L;
        else
        fprintf('Invalid method. Choose 1 for eig or 2 for chol.\n');
end
        Z = res*Y;
        X = Z + m.*ones(n,N);
    fprintf('\tThe size of matrix X is %dx%d.\n', n, N);
end