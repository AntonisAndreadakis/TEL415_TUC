%TEL415 - ERGASIA 3
%ARI8MOS OMADAS: LAB41544983
%Matsatsos Ioannis 2013030148
%Andreadakis Antonis 2013030059
clc; clear; close all;


%% Askhsh 1
%% erwthma A
T = 0.5;    % secs
N = 18;     % bits {-1, +1}

Ts = 10^-5;     % periodos deigmatolh4ias
Fs = 1/Ts;      % sixnothta deigmatolh4ias

t0 = [0:Ts:T-Ts];       % gia 1 periodo
t = [0:Ts:(N+2)*T-Ts];  % gia olo to shma

A = 9.798;  % Aplitude for T=0.5

% creating pulse p
pa = A * t0(t0<T/2);
pb = A * (T-t0(t0>=T/2));
p = [pa pb];

% creating signal x baseband
x = [];
for i = 1:18
    s = 2 * (rand > 0.5) - 1;
    p_new = sqrt(10) * s * p;
    x = [x p_new];
end
x = [x zeros(1, 2 * length(p))];    % last two periods are zero

X = fftshift(fft(x)) * Ts;          % Fourier of x
K = length(x);
F = [-K/2:K/2-1] * Fs/K;

%figure;
%plot(t, x);
figure
plot(F, abs(X))
xlim([-100 100]); xlabel('F'); legend('|X(F)|')

%% erwthma B
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
f1=find(filter_lowpass);
X_filt = X .* filter_lowpass;   % X filtered

x_filt = ifft(fftshift(X_filt)) * Fs;       % Inverse Fourier of X_filt

figure
plot(t, x, t, x_filt,'--')
xlabel('t'); legend('x(t)','x_f_i_l_t(t)')

% x bandpass
Fc = 1000;          % Hz - carrier frequency
x_BP = x_filt .* cos(2 * pi * Fc * t);      % modulation with 1000 Hz carrier

figure;
plot(t, x_BP)
xlabel('t'); legend('x_B_P(t)')

X_BP = fftshift(fft(x_BP)) * Ts;        % Fourier of x bandpass

figure
plot(F, abs(X_BP))
xlim([900 1100]); xlabel('F'); legend('|X_B_P(F)|')

%% erwthma C
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

%% erwthma D
% move positive frequencies' content of Y_BP_filt to [-B,+B]
%Y = [zeros(1, 499900) Y_BP_filt(509901:510100) zeros(1, 499900)];
Y = [zeros(1, K/2-B_points) Y_BP_filt(K/2-B_points+B_points^2+1:K/2+B_points^2+B_points) zeros(1, K/2-B_points)];

y = ifft(fftshift(Y)) * Fs;     % Inverse Fourier of Y

h = sum(1/2 * al .* exp(-1i*2*pi*Fc.*tl));

figure
plot(t, x_filt, t, real(y./h))
hold on

% erwthma (C, D), for-loop for all delays
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

%% erwthma E
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