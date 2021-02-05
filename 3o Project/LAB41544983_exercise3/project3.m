%TEL415 - ERGASIA 3
%ARI8MOS OMADAS: LAB41544983
%Matsatsos Ioannis 2013030148
%Andreadakis Antonis 2013030059
        %    ASKHSH 1  %
        clear; clc; close all;
        fprintf('\n######### Askhsh 1: #########\n');
% EROTIMA A:
N = 18;  % bits {-1, +1}
T = 0.5; % secs
E = 10;
Ts = 10^-5; % periodos deigmatolhpsias
Fs = 1/Ts;  % sixnothta deigmatolhpsias
t0 = [0:Ts:T-Ts];      % gia 1 periodo
t = [0:Ts:(N+2)*T-Ts]; % gia olo to shma
Ep = 1; % energeia palmoy. Omws Ep = olokliroma (p^2)(t)dt apo
        % mhden ews T. Stin periptwsi mas exoyme T = 0.5 sec
        % Lunontas to parapanw oloklirwma exoyme 0.014167*(Á^2) = 1, diladi
        % A^2 = 95.999628003036854, opote telika A = 9.798 (afou A>0).
A = 9.798;  % Aplitude for T=0.5
% creating pulse p
pa = A * t0(t0<T/2);
pb = A * (T-t0(t0>=T/2));
p = [pa pb];
% creating signal x baseband
x = [];
for i = 1:N
    s = 2 * (rand > 0.5) - 1;
    p_new = sqrt(E)*s*p;
    x = [x p_new];
end
x = [x zeros(1, 2 * length(p))];    % last two periods are zero
X = fftshift(fft(x)) * Ts;          % Fourier of x
K = length(x);
F = [-K/2:K/2-1] * Fs/K;
figure;
plot(t, x);

figure
plot(F, abs(X))
xlim([-100 100]); xlabel('F'); legend('|X(F)|')

% EROTIMA B
B = 10;    % Hz
B_points = round(B*K*Ts);
% 1os tropos - slow
%for i = 1:length(F)
%    if abs(F(i)) <= B
%        X_filt(i) = X(i);
%    else
%        X_filt(i) = 0;
%    end
%end

% 2os tropos - fast
filter_lowpass = [zeros(1, K/2-B_points) ones(1, 2*B_points) zeros(1, K/2-B_points)];
%f1=find(filter_lowpass);
X_filt = X .* filter_lowpass;   % X filtered
x_filt = ifft(fftshift(X_filt)) * Fs; % Inverse Fourier of X_filt

figure
plot(t, x, t, x_filt,'--')
xlabel('t'); legend('x(t)','x_f_i_l_t(t)')

% x bandpass
Fc = 1000;          % Hz - carrier frequency
x_BP = x_filt .* cos(2 * pi * Fc * t); % modulation with 1000 Hz carrier

figure;
plot(t, x_BP)
xlabel('t'); legend('x_B_P(t)')

X_BP = fftshift(fft(x_BP)) * Ts;  % Fourier of x bandpass

figure
plot(F, abs(X_BP))
xlim([900 1100]); xlabel('F'); legend('|X_B_P(F)|')

% EROTIMA C
L = 100;        % number of arrays
al = rand(1, L);
max_delay = 0.1 * T;
tl = max_delay * rand(1, L);
% create y bandpass
for i = 1:L
    tl_to_zeros = round(tl(i) / Ts);
    x_BP_delayed_tl = [zeros(1,tl_to_zeros) x_BP(1:end-tl_to_zeros)];
    y_temp(i, :) = al(i) .* x_BP_delayed_tl;
end
y_BP = sum(y_temp);
Y_BP = fftshift(fft(y_BP)) * Ts;        % Fourier of y bandpass
% 1os tropos - slow
%for i = 1:length(F)
%    if abs(F(i)-Fc) <= B
%        Y_BP_filt(i) = Y_BP(i);
%    elseif abs(F(i)+Fc) <= B
%        Y_BP_filt(i) = Y_BP(i);
%    else
%        Y_BP_filt(i) = 0;
%    end
%end

% 2os tropos - fast
%filter_bandpass = [zeros(1,489900) ones(1,200) zeros(1,19800) ones(1,200) zeros(1,489900)];
filter_bandpass = [zeros(1, K/2-B_points^2-B_points) ones(1, 2*B_points) zeros(1, 2*B_points^2-2*B_points) ones(1, 2*B_points) zeros(1, K/2-B_points^2-B_points)];
Y_BP_filt = Y_BP .* filter_bandpass;        % bandpass filter Y
y_BP_filt = ifft(fftshift(Y_BP_filt)) * Fs;     % inverse Fourier of Y bandpass

figure;
plot(F, abs(Y_BP))
xlim([900 1100]); xlabel('F'); legend('|Y_B_P(F)|')

figure
plot(F, abs(Y_BP_filt))
xlim([900 1100]); xlabel('F'); legend('|Y_B_P_,_f_i_l_t(F)|')

% EROTIMA D
% move positive frequencies' content of Y_BP_filt to [-B,+B]
%Y = [zeros(1, 499900) Y_BP_filt(509901:510100) zeros(1, 499900)];
Y = [zeros(1, K/2-B_points) Y_BP_filt(K/2-B_points+B_points^2+1:K/2+B_points^2+B_points) zeros(1, K/2-B_points)];
y = ifft(fftshift(Y)) * Fs;     % Inverse Fourier of Y

h = sum(1/2 * al .* exp(-1i*2*pi*Fc.*tl));

figure
plot(t, x_filt, t, real(y./h))
hold on

% EROTIMA (C, D), for-loop for all delays
for delay = [0.2 0.3 0.5 1.5]
    max_delay = delay * T;
    tl = max_delay * rand(1, L);
    for i = 1:L
        tl_to_zeros = round(tl(i) / Ts);
        x_BP_delayed_tl = [zeros(1,tl_to_zeros) x_BP(1:end-tl_to_zeros)];
        y_temp(i, :) = al(i) .* x_BP_delayed_tl;
    end
    y_BP = sum(y_temp);
    Y_BP = fftshift(fft(y_BP)) * Ts;
    filter_bandpass = [zeros(1, K/2-B_points^2-B_points) ones(1, 2*B_points) zeros(1, 2*B_points^2-2*B_points) ones(1, 2*B_points) zeros(1, K/2-B_points^2-B_points)];
    Y_BP_filt = Y_BP .* filter_bandpass;
    y_BP_filt = ifft(fftshift(Y_BP_filt)) * Fs;
    Y = [zeros(1, K/2-B_points) Y_BP_filt(K/2-B_points+B_points^2+1:K/2+B_points^2+B_points) zeros(1, K/2-B_points)];
    y = ifft(fftshift(Y)) * Fs;
    h = sum(1/2 * al .* exp(-1i*2*pi*Fc.*tl));
    plot(t, real(y./h));
end
xlabel('t'); legend('x_f_i_l_t','y, delay=0.1*T','y, delay=0.2*T','y, delay=0.3*T','y, delay=0.5*T','y, delay=1.5*T')

% EROTIMA E
L = 400;
for i = 1 : 10^6
    al = rand(1, L);
    max_delay = 0.1 * T;
    tl = max_delay * rand(1, L);
    h(i) = sum(1/2 * al .* exp(-1i*2*pi*Fc.*tl));
end

hr = real(h); hi = imag(h);

figure
histogram2(hr,hi);
xlabel('real(h)'); ylabel('imag(h)'); title('pdf of h')

% repeat above for max_delay = 0
for i = 1 : 10^6
    al = rand(1, L);
    max_delay = 0 * T;
    tl = max_delay * rand(1, L);
    h0(i) = sum(1/2 * al .* exp(-1i*2*pi*Fc.*tl));
end

h0r = real(h0); h0i = imag(h0);

figure
histogram2(h0r,h0i);
xlabel('real(h)'); ylabel('imag(h)'); title('pdf of h, max-delay = 0')
pause;
            %    ASKHSH 2  %
            clear; clc; close all;
            fprintf('\n######### Askhsh 2: #########\n');
% EROTIMA A:
load('real_data.mat')
fprintf('\n');
fprintf('\tWe will use a matrix decomposition method for real data.\n');
fprintf('\tWhich method would you like to use?\n');
fprintf('\tPress 1 for eigendecomposition or 2 for Cholesky:\n');
method = input(' ');
N = 1000; % N is number of independent vectors
X = real_colored_noise(m, C, N, method);
real_mean = mean(X,2);
fprintf('\tmean vector:\n');
disp(real_mean);
real_C = cov(X');
fprintf('\tcovariance matrix:\n');
disp(real_C);
                   
% EROTIMA B:
load('complex_data.mat')
fprintf('\n');
fprintf('\tWe will use a matrix decomposition method for complex data.\n');
fprintf('\tWhich method would you like to use?\n');
fprintf('\tPress 1 for eigendecomposition or 2 for Cholesky:\n');
method = input(' ');
N = 1000; % N is number of independent vectors
X2 = complex_colored_noise(m,C,N,method);
complex_mean = mean(X2,2);
fprintf('\tmean vector:\n');
disp(complex_mean);
complex_C = cov(X2');
fprintf('\tcovariance matrix:\n');
disp(complex_C);
pause;

          %    ASKHSH 3  %
          clear; clc; close all;
          fprintf('\n######### Askhsh 3: #########\n');
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
pause;

                %    ASKHSH 4  %
        clear; clc; close all;
        fprintf('\n######### Askhsh 4: #########\n');
        % we want DoA1...DoA5 as BPSK
m = 20; % stoixeia gia to antenna array
n = 5; % 5 piges
DoA = [-75 -60 0 60 75]; % oi antistoixes gwnies me to array
r = length(DoA); % Total number of signals

% H apostasi didoxikwn stoixeiwn tou array: d = lambda/2,
% lambda: mikos kymatos, ara d/lambda = 1/2. Opote anti gia
% d/lambda bazw 1/2:

% EROTIMA 1 A:
% given Power in dB:
P_dB1 = [7 5 7 5 7];
P1 = 10.^(P_dB1/10); % convert to watts
 % create a matrix with m rows and 1 column
for i = 0 : (m-1)
    % create a(theta):
    a1(i+1,1) = exp(-1i*pi*i*sin(deg2rad(-75)));
    a2(i+1,1) = exp(-1i*pi*i*sin(deg2rad(-60)));
    a3(i+1,1) = exp(-1i*pi*i*sin(deg2rad(0)));
    a4(i+1,1) = exp(-1i*pi*i*sin(deg2rad(60)));
    a5(i+1,1) = exp(-1i*pi*i*sin(deg2rad(75)));
end
% no need to compute noise..
% Estimating the covariance matrix of the sensor array %%%%%%%
R = (P1(1)*(a1*a1')+P1(2)*(a2*a2')+P1(3)*(a3*a3')+P1(4)*(a4*a4')+P1(5)*(a5*a5'))+ eye(m);      
[U,S,V] = svd(R); %Compute eigendecomposition of covariance matrix
Qn = U(:,r+1:m); %Get the noise eigenvectors
IR = inv(R); %Inverse of covariance matrix
theta = -90:1:90;
% Compute steering vectors corresponding values in angles:
for i = 1:length(theta)
    % dianismata a, gia ton ypologismo twn ektimitwn:
    for k = 0 : m-1
        % to 2 me to 1/2 aplopoieitai opws eidame parapanw
        % afoy d/lambda = 1/2 kai 2*pi*d/lambda = pi
        angles(k+1,1) = exp(-1i*(pi)*k*sin(deg2rad(theta(i))));
    end 
    %MF algorithm
    mf(i) = real(angles'*R*angles);
    %MVDR algorithm
    mvdr(i)= real(1/(angles'*IR*angles));
    % MUSIC algorithm
    music_spectrum(i) = 1/(norm(Qn'*angles).^2);
end
mf = mf./max(mf);
mvdr = mvdr./max(mvdr);
music_spectrum = music_spectrum./max(music_spectrum);
    % ESPRIT algorithm:    
    var = diag(S); % take all diag elements of S
    var = mean(var(n+1:m)); % chose the small values
    Rr = (R-var*eye(m)); % minus (sigma^2)*I as given in page 66 slide19
    Resp = [Rr(1:m-1,:); Rr(2:m,:)]; % create matrix for esprit
    [Ue, Se, Ve] = svd(Resp); % use decomposition on that matrix
    U1 = Ue(1:38/2,1:n); % separate Ue equaly in 2 matrices
    U2 = Ue(38/2+1:38,1:n);
    R1 = U1'*U1;
    R2 = U1'*U2;
    [E, D] = eig(inv(R1)*R2); % again use decomposition on new matix
    e = diag(D); % now we have values corresponding to our angles
    omega = -1i*log((e));%take natural log to find omega from exp(j*omega) as given in page 66 slide19
    theta_esprit = real(asind(-omega/pi)); % at this point we have our angles
figure(1)
plot(theta,mf,'b')
xlabel('\theta')
ylabel(' P(\theta) ')
hold on
plot(theta,mvdr,'r')
plot(theta,music_spectrum,'m')
% conclude DoA1...DoA5 using ESPRIT:
% normalize every power spec from 0 to 1:
% draw one line per column for each value of theta_esprit:
line([theta_esprit(1) theta_esprit(1)],[0 1],'Color','k');
line([theta_esprit(2) theta_esprit(2)],[0 1],'Color','k');
line([theta_esprit(3) theta_esprit(3)],[0 1],'Color','k');
line([theta_esprit(4) theta_esprit(4)],[0 1],'Color','k');
line([theta_esprit(5) theta_esprit(5)],[0 1],'Color','k');
% draw true values of angles:
% normalize every power spec from 0 to 1:
plot([75 75],[0 1],'--g');
plot([-75 -75],[0 1],'--g');
plot([60 60],[0 1],'--g');
plot([-60 -60],[0 1],'--g');
plot([0 0],[0 1],'--g');
title('DoA Estimation based on MF, MVDR and MUSIC Algorithms')
legend('MF','MVDR','MUSIC','ESPRIT', 'true');
grid on
axis([-90 90 0 1])

%Erotima B,C,D (change only variable N from 50 to 200 and 1000):
A = [a1 a2 a3 a4 a5];
P2 = diag(sqrt(P1));
N = 1000;
b = sign(randn(n,N)); % transmitted bit
% white Complex Gaussian noise:
w = (randn(m, N)+(1i*randn(m, N)))/sqrt(2);
Y = A*P2*b+w; % our signal as given at equation
R2 = (1/N).*(Y*Y'); %correlation matrix
%decomposition:
[U,S,V] = svd(R2);
Qn = U(:,n+1:m);
theta = -90:1:90;
for i = 1 : length(theta)
    for k = 0 : m-1
        angles(k+1,1) = exp(-1i*(pi)*k*sin(deg2rad(theta(i))));
    end
    %MF algorithm
    mf(i) = real(angles'*R*angles);
    %MVDR algorithm
    mvdr(i)= real(1/(angles'*IR*angles));
    % MUSIC algorithm
    music_spectrum(i) = 1/(norm(Qn'*angles).^2);
end
mf = mf./max(mf);
mvdr = mvdr./max(mvdr);
music_spectrum = music_spectrum./max(music_spectrum);
    % ESPRIT algorithm:    
    var = diag(S); % take all diag elements of S
    var = mean(var(n+1:m)); % chose the small values
    Rr = (R-var*eye(m)); % minus (sigma^2)*I as given in page 66 slide19
    Resp = [Rr(1:m-1,:); Rr(2:m,:)]; % create matrix for esprit
    [Ue, Se, Ve] = svd(Resp); % use decomposition on that matrix
    U1 = Ue(1:38/2,1:n); % separate Ue equaly in 2 matrices
    U2 = Ue(38/2+1:38,1:n);
    R1 = U1'*U1;
    R2 = U1'*U2;
    [E, D] = eig(inv(R1)*R2); % again use decomposition on new matix
    e = diag(D); % now we have values corresponding to our angles
    omega = -1i*log((e)); %take natural log to find omega from exp(j*omega) as given in page 66 slide19
    theta_esprit = real(asind(-omega/pi)); % at this point we have our angles
figure(2)
plot(theta,mf,'b')
xlabel('\theta')
ylabel(' P(\theta) ')
hold on
plot(theta,mvdr,'r')
plot(theta,music_spectrum,'m')
% conclude DoA1...DoA5 using ESPRIT:
% normalize every power spec from 0 to 1:
% draw one line per column for each value of theta_esprit:
line([theta_esprit(1) theta_esprit(1)],[0 1],'Color','k');
line([theta_esprit(2) theta_esprit(2)],[0 1],'Color','k');
line([theta_esprit(3) theta_esprit(3)],[0 1],'Color','k');
line([theta_esprit(4) theta_esprit(4)],[0 1],'Color','k');
line([theta_esprit(5) theta_esprit(5)],[0 1],'Color','k');
% draw true values of angles:
% normalize every power spec from 0 to 1:
plot([75 75],[0 1],'--g');
plot([-75 -75],[0 1],'--g');
plot([60 60],[0 1],'--g');
plot([-60 -60],[0 1],'--g');
plot([0 0],[0 1],'--g');
title('P = [7 5 7 5 7] N = 1000')
legend('MF','MVDR','MUSIC','ESPRIT', 'true');
grid on
axis([-90 90 0 1])
%pause;

% Erotima 2 A:
% given Power in dB:
P_dB2 = [-7 -5 -7 -5 -7];
P3 = 10.^(P_dB2/10); % convert to watts
 % create a matrix with m rows and 1 column
for i = 0 : (m-1)
    % create a(theta):
    a1(i+1,1) = exp(-1i*pi*i*sin(deg2rad(-75)));
    a2(i+1,1) = exp(-1i*pi*i*sin(deg2rad(-60)));
    a3(i+1,1) = exp(-1i*pi*i*sin(deg2rad(0)));
    a4(i+1,1) = exp(-1i*pi*i*sin(deg2rad(60)));
    a5(i+1,1) = exp(-1i*pi*i*sin(deg2rad(75)));
end
% no need to compute noise..
% Estimating the covariance matrix of the sensor array %%%%%%%
R = ((P3(1))*(a1*a1')+(P3(2))*(a2*a2')+(P3(3))*(a3*a3')+(P3(4))*(a4*a4')+(P3(5))*(a5*a5'))+ eye(m);      
[U,S,V] = svd(R); %Compute eigendecomposition of covariance matrix
Qn = U(:,r+1:m); %Get the noise eigenvectors
IR = inv(R); %Inverse of covariance matrix
theta = -90:1:90;
% Compute steering vectors corresponding values in angles:
for i = 1:length(theta)
    % dianismata a, gia ton ypologismo twn ektimitwn:
    for k = 0 : m-1
        % to 2 me to 1/2 aplopoieitai opws eidame parapanw
        % afoy d/lambda = 1/2 kai 2*pi*d/lambda = pi
        angles(k+1,1) = exp(-1i*(pi)*k*sin(deg2rad(theta(i))));
    end 
    %MF algorithm
    mf(i) = real(angles'*R*angles);
    %MVDR algorithm
    mvdr(i)= real(1/(angles'*IR*angles));
    % MUSIC algorithm
    music_spectrum(i) = 1/(norm(Qn'*angles).^2);
end
mf = mf./max(mf);
mvdr = mvdr./max(mvdr);
music_spectrum = music_spectrum./max(music_spectrum);
    % ESPRIT algorithm:    
    var = diag(S); % take all diag elements of S
    var = mean(var(n+1:m)); % chose the small values
    Rr = (R-var*eye(m));  % minus (sigma^2)*I as given in page 66 slide19
    Resp = [Rr(1:m-1,:); Rr(2:m,:)]; % create matrix for esprit
    [Ue, Se, Ve] = svd(Resp); % use decomposition on that matrix
    U1 = Ue(1:38/2,1:n); % separate Ue equaly in 2 matrices
    U2 = Ue(38/2+1:38,1:n); 
    R1 = U1'*U1; 
    R2 = U1'*U2; 
    [E, D] = eig(inv(R1)*R2); % again use decomposition on new matix
    e = diag(D); % now we have values corresponding to our angles
    omega = -1i*log((e)); %take natural log to find omega from exp(j*omega) as given in page 66 slide19
    theta_esprit = real(asind(-omega/pi)); % at this point we have our angles
figure(3)
plot(theta,mf,'b')
xlabel('\theta')
ylabel(' P(\theta) ')
hold on
plot(theta,mvdr,'r')
plot(theta,music_spectrum,'m')
% conclude DoA1...DoA5 using ESPRIT:
% normalize every power spec from 0 to 1:
% draw one line per column for each value of theta_esprit:
line([theta_esprit(1) theta_esprit(1)],[0 1],'Color','k');
line([theta_esprit(2) theta_esprit(2)],[0 1],'Color','k');
line([theta_esprit(3) theta_esprit(3)],[0 1],'Color','k');
line([theta_esprit(4) theta_esprit(4)],[0 1],'Color','k');
line([theta_esprit(5) theta_esprit(5)],[0 1],'Color','k');
% draw true values of angles:
% normalize every power spec from 0 to 1:
plot([75 75],[0 1],'--g');
plot([-75 -75],[0 1],'--g');
plot([60 60],[0 1],'--g');
plot([-60 -60],[0 1],'--g');
plot([0 0],[0 1],'--g');
title('DoA Estimation based on MF, MVDR and MUSIC Algorithms')
legend('MF','MVDR','MUSIC','ESPRIT', 'true');
grid on
axis([-90 90 0 1])

%Erotima B,C,D (change only variable N from 50 to 200 and 1000):
A = [a1 a2 a3 a4 a5];
P2 = diag(sqrt(P3));
N = 1000;
b = sign(randn(n,N));
% white Complex Gaussian noise:
w = (randn(m, N)+(1i*randn(m, N)))/sqrt(2);
Y = A*P2*b+w;
R2 = (1/N).*(Y*Y');
%decomposition:
[U,S,V] = svd(R2);
Qn = U(:,n+1:m);
theta = -90:1:90;
for i = 1 : length(theta)
    for k = 0 : m-1
        angles(k+1,1) = exp(-1i*(pi)*k*sin(deg2rad(theta(i))));
    end
    %MF algorithm
    mf(i) = real(angles'*R*angles);
    %MVDR algorithm
    mvdr(i)= real(1/(angles'*IR*angles));
    % MUSIC algorithm
    music_spectrum(i) = 1/(norm(Qn'*angles).^2);
end
mf = mf./max(mf);
mvdr = mvdr./max(mvdr);
music_spectrum = music_spectrum./max(music_spectrum);
    % ESPRIT algorithm:    
    var = diag(S); % take all diag elements of S
    var = mean(var(n+1:m)); % chose the small values
    Rr = (R-var*eye(m)); % minus (sigma^2)*I as given in page 66 slide19
    Resp = [Rr(1:m-1,:); Rr(2:m,:)]; % create matrix for esprit
    [Ue, Se, Ve] = svd(Resp); % use decomposition on that matrix
    U1 = Ue(1:38/2,1:n); % separate Ue equaly in 2 matrices
    U2 = Ue(38/2+1:38,1:n);
    R1 = U1'*U1;
    R2 = U1'*U2;
    [E, D] = eig(inv(R1)*R2); % again use decomposition on new matix
    e = diag(D); % now we have values corresponding to our angles
    omega = -1i*log((e));%take natural log to find omega from exp(j*omega) as given in page 66 slide19
    theta_esprit = real(asind(-omega/pi)); % at this point we have our angles
figure(4)
plot(theta,mf,'b')
xlabel('\theta')
ylabel(' P(\theta) ')
hold on
plot(theta,mvdr,'r')
plot(theta,music_spectrum,'m')
% conclude DoA1...DoA5 using ESPRIT:
% normalize every power spec from 0 to 1:
% draw one line per column for each value of theta_esprit:
line([theta_esprit(1) theta_esprit(1)],[0 1],'Color','k'); 
line([theta_esprit(2) theta_esprit(2)],[0 1],'Color','k');
line([theta_esprit(3) theta_esprit(3)],[0 1],'Color','k');
line([theta_esprit(4) theta_esprit(4)],[0 1],'Color','k');
line([theta_esprit(5) theta_esprit(5)],[0 1],'Color','k');
% draw true values of angles:
% normalize every power spec from 0 to 1:
plot([75 75],[0 1],'--g');
plot([-75 -75],[0 1],'--g');
plot([60 60],[0 1],'--g');
plot([-60 -60],[0 1],'--g');
plot([0 0],[0 1],'--g');
title('P = [-7 -5 -7 -5 -7] N = 1000')
legend('MF','MVDR','MUSIC','ESPRIT', 'true');
grid on
axis([-90 90 0 1])
%pause;





% FUNCTIONS FOR 2nd EXERCISE:

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