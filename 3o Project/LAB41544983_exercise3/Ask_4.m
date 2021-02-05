%TEL415 - ERGASIA 3
%ARI8MOS OMADAS: LAB41544983
%Matsatsos Ioannis 2013030148
%Andreadakis Antonis 2013030059
clc; clear all; close all;

%% ASKHSH 4
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
    omega = -1i*log((e)); %take natural log to find omega from exp(j*omega) as given in page 66 slide19
    theta_esprit = real(asind(-omega/pi)); % at this point we have our angles
figure(1)
plot(theta,mf,'b')
xlabel('\theta')
ylabel(' P(\theta) ')
hold on
plot(theta,mvdr,'r')
plot(theta,music_spectrum,'m')
% normalize every power spec from 0 to 1.
line([theta_esprit(1) theta_esprit(1)],[0 1],'Color','k');
line([theta_esprit(2) theta_esprit(2)],[0 1],'Color','k');
line([theta_esprit(3) theta_esprit(3)],[0 1],'Color','k');
line([theta_esprit(4) theta_esprit(4)],[0 1],'Color','k');
line([theta_esprit(5) theta_esprit(5)],[0 1],'Color','k');
% conclude DoA1...DoA5 using ESPRIT:
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
% normalize every power spec from 0 to 1.
line([theta_esprit(1) theta_esprit(1)],[0 1],'Color','k');
line([theta_esprit(2) theta_esprit(2)],[0 1],'Color','k');
line([theta_esprit(3) theta_esprit(3)],[0 1],'Color','k');
line([theta_esprit(4) theta_esprit(4)],[0 1],'Color','k');
line([theta_esprit(5) theta_esprit(5)],[0 1],'Color','k');
% conclude DoA1...DoA5 using ESPRIT:
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
% normalize every power spec from 0 to 1.
line([theta_esprit(1) theta_esprit(1)],[0 1],'Color','k');
line([theta_esprit(2) theta_esprit(2)],[0 1],'Color','k');
line([theta_esprit(3) theta_esprit(3)],[0 1],'Color','k');
line([theta_esprit(4) theta_esprit(4)],[0 1],'Color','k');
line([theta_esprit(5) theta_esprit(5)],[0 1],'Color','k');
% conclude DoA1...DoA5 using ESPRIT:
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
    omega = -1i*log((e)); %take natural log to find omega from exp(j*omega) as given in page 66 slide19
    theta_esprit = real(asind(-omega/pi)); % at this point we have our angles
figure(4)
plot(theta,mf,'b')
xlabel('\theta')
ylabel(' P(\theta) ')
hold on
plot(theta,mvdr,'r')
plot(theta,music_spectrum,'m')
% normalize every power spec from 0 to 1.
line([theta_esprit(1) theta_esprit(1)],[0 1],'Color','k');
line([theta_esprit(2) theta_esprit(2)],[0 1],'Color','k');
line([theta_esprit(3) theta_esprit(3)],[0 1],'Color','k');
line([theta_esprit(4) theta_esprit(4)],[0 1],'Color','k');
line([theta_esprit(5) theta_esprit(5)],[0 1],'Color','k');
% conclude DoA1...DoA5 using ESPRIT:
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
