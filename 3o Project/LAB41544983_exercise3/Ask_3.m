%TEL415 - ERGASIA 3
%ARI8MOS OMADAS: LAB41544983
%Matsatsos Ioannis 2013030148
%Andreadakis Antonis 2013030059
clc; clear all; close all;

%% ASKHSH 3
m = 31;         % bit-signature length
n = 10;         % number of users
N = 10^6;       % bits total per user
K = 200;        % number of samples for SMI, SA-SMI filters

%s = (2 * (rand(m, n) > 0.5) - 1)/sqrt(m);  % signatures of users(manual)
load('signatures.mat');
s=S;
P_dB = [0 5 5 5 10 10 10 15 15 15];        % power of all users in dB

b = 2 * (rand(n,N) > 0.5) - 1;       % bits of n users for N transmissions

noise = randn(m, N);        % white gaussian noise with P=0dB;

% erwthmata 1,2,3
for P1 = 0 : 20
    P_dB(1) = P1;       % power of user 1 in dB
    P = 10.^(P_dB/10);      % power in power units
    
    % received signal by base station
    x = s * diag(sqrt(P)) * b + noise;
    
    % FILTERS
    % MF filter
    w_mf = s(:,1);
    % MMSE filter
    Rx = s * diag(P) * s' + eye(m);
    w_mmse = inv(Rx) * s(:,1);
    % ZF filter
    w_zf = s * inv(s' * s);
    w_zf = w_zf(:,1);
    % SMI filter
    Rx_smi = (1/K).*(x(:,1:K) * x(:,1:K)');
    w_smi = inv(Rx_smi) * s(:,1);
    % SA-SMI filter
    z = x(:,1:K) - (s(:,1) * diag(sqrt(P(1))) * b(1,1:K));
    Rz_sasmi = (1/K).*(z * z');
    w_sasmi = inv(Rz_sasmi) * s(:,1);
    % SA-SMI loop filter
    y_loop = w_smi' * x;    % use smi filter first
    b_hat_loop = sign(y_loop);      % decide
    while true
        z_hat = x(:,1:K) - (s(:,1) * diag(sqrt(P(1))) * b_hat_loop(:,1:K));   % estimate z
        Rz_hat = (1/K).*(z_hat * z_hat');   % create Rz based on z est
        w_sasmi_loop = inv(Rz_hat) * s(:,1);    % create sa_smi filter
        y_loop_new = w_sasmi_loop' * x;    % filter
        b_hat_loop_new = sign(y_loop_new);      % decide
    
        if (b_hat_loop == b_hat_loop_new)       % if two btb decisions are the same, break inf loop
            break;
        else
            b_hat_loop = b_hat_loop_new;        % else old decision = new decision
        end
    end
    
    % FILTERING
    y_mf = w_mf' * x;
    y_mmse = w_mmse' * x;
    y_zf = w_zf' * x;
    y_smi = w_smi' * x;
    y_sasmi = w_sasmi' * x;
    
    % DECIDE
    b_hat_mf = sign(y_mf);
    b_hat_mmse = sign(y_mmse);
    b_hat_zf = sign(y_zf);
    b_hat_smi = sign(y_smi);
    b_hat_sasmi = sign(y_sasmi);
    
    % BER calculation
    BER_mf(P1+1) = size(find([(b(1,:) - b_hat_mf)]),2)/N;
    BER_mmse(P1+1) = size(find([(b(1,:) - b_hat_mmse)]),2)/N;
    BER_zf(P1+1) = size(find([(b(1,:) - b_hat_zf)]),2)/N;
    BER_smi(P1+1) = size(find([(b(1,:) - b_hat_smi)]),2)/N;
    BER_sasmi(P1+1) = size(find([(b(1,:) - b_hat_sasmi)]),2)/N;
    BER_sasmi_loop(P1+1) = size(find([(b(1,:) - b_hat_loop)]),2)/N;
end

figure
semilogy([0:20], BER_mf, 'b')
hold on;
semilogy([0:20], BER_mmse, 'r')
hold on;
semilogy([0:20], BER_zf, 'black')
hold on;
semilogy([0:20], BER_smi, 'green')
hold on;
semilogy([0:20], BER_sasmi, 'm')
hold on;
semilogy([0:20], BER_sasmi_loop, 'c')

% erwthma 4
for P1 = 0 : 20
    P_dB(1) = P1;       % power of user 1 in dB
    P = 10.^(P_dB/10);      % power in power units
    
    % optimising user's 1 signature (real number)
    H = (s * diag(sqrt(P)) * s')./norm(s*s') + eye(m);     % channel matrix
    Rz = s(:,2:end) * diag(P(1,2:end)) * s(:,2:end)' + eye(m);
    [Q,L] = eig(H' * inv(Rz) * H);
    [maxL, maxidx] = max(sum(L));
    %norm(Q(:,maxidx))
    s(:,1) = Q(:,maxidx);
    
    % received signal by base station
    x = s * diag(sqrt(P)) * b + noise;
    
    % FILTERS
    % MF filter
    w_mf = s(:,1);
    % MMSE filter
    Rx = s * diag(P) * s' + eye(m);
    w_mmse = inv(Rx) * s(:,1);
    % ZF filter
    w_zf = s * inv(s' * s);
    w_zf = w_zf(:,1);
    % SMI filter
    Rx_smi = (1/K).*(x(:,1:K) * x(:,1:K)');
    w_smi = inv(Rx_smi) * s(:,1);
    % SA-SMI filter
    z = x(:,1:K) - (s(:,1) * diag(sqrt(P(1))) * b(1,1:K));
    Rz_sasmi = (1/K).*(z * z');
    w_sasmi = inv(Rz_sasmi) * s(:,1);
    
    % FILTERING
    y_mf = w_mf' * x;
    y_mmse = w_mmse' * x;
    y_zf = w_zf' * x;
    y_smi = w_smi' * x;
    y_sasmi = w_sasmi' * x;
    
    % DECIDE
    b_hat_mf = sign(y_mf);
    b_hat_mmse = sign(y_mmse);
    b_hat_zf = sign(y_zf);
    b_hat_smi = sign(y_smi);
    b_hat_sasmi = sign(y_sasmi);
    
    % BER calculation
    BER_mf(P1+1) = size(find([(b(1,:) - b_hat_mf)]),2)/N;
    BER_mmse(P1+1) = size(find([(b(1,:) - b_hat_mmse)]),2)/N;
    BER_zf(P1+1) = size(find([(b(1,:) - b_hat_zf)]),2)/N;
    BER_smi(P1+1) = size(find([(b(1,:) - b_hat_smi)]),2)/N;
    BER_sasmi(P1+1) = size(find([(b(1,:) - b_hat_sasmi)]),2)/N;
end

semilogy([0:20], BER_mf, 'b-*')
hold on;
semilogy([0:20], BER_mmse, 'r-o')
hold on;
semilogy([0:20], BER_zf, 'black--')
hold on;
semilogy([0:20], BER_smi, 'green--')
hold on;
semilogy([0:20], BER_sasmi, 'm--')
xlabel('SNR(dB)'); ylabel('BER'); legend('MF','MMSE','ZF','SMI','SA-SMI','SA-SMI-loop','MF(opt. s_1)','MMSE(opt. s_1)','ZF(opt. s_1)','SMI(opt. s_1)','SA-SMI(opt. s_1)');