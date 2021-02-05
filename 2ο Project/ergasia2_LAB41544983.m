%TEL415 - ERGASIA 2
%ARI8MOS OMADAS: LAB41544983
%Matsatsos Ioannis 2013030148
%Andreadakis Antonis 2013030059
clc; clear; close all;


%% Askhsh 1
fprintf('\n######### Askhsh 1 #########\n');
%test myQR(A) with a random Imaginary 4x6 matrix
%m<n=rank
m=4; n=6;
A=rand(m,n)+1i*rand(m,n);

% 2 checks
[Q,R]=myQR(A);
fprintf('\t Max of max of abs(A-Q*R): \n');
myresult=max(max(abs(A-Q*R)));
display(myresult);

fprintf('\t Max of max of abs(Q"*Q-I): \n');
I=eye(size(A,1),rank(A));
max(max(abs(Q'*Q-I)))

% compare with matlab function
fprintf('\t Difference of myQR to matlab function: \n');
[Qm,Rm]=qr(A);
matlabresult=max(max(abs(A-Qm*Rm)));
difference=myresult-matlabresult;
display(difference);

fprintf('\nPaused. Press enter to continue.\n');
pause;

%% Askhsh 3
clc; close all; clear all;
fprintf('\n######### Askhsh 3 #########\n');

n=10000;
X = randn(n,1);
Y = randn(n,1);
for i=1:n
    W(i,1)=(2*(randi(2)>1)-1)*X(i);
end
%%%%%%%% Checking if R(Z) and I(Z) are uncorrelated %%%%%%%
a1=corr(X,X);
%a2=cov(X,X);
if (round(a1)==1 || round(a1)==-1)
    fprintf('The R(Z)=X and I(Z)=X are correlated with r = %.0f\n',a1);
elseif (round(a1)==0)
    fprintf('The R(Z)=X and I(Z)=X are uncorrelated\n');
end
b1=corr(X,Y);
%b2=cov(X,Y);
if (round(b1)==1 || round(b1)==-1)
    fprintf('The R(Z)=X and I(Z)=Y are correlated with r = %d\n',b1);
elseif (round(b1)==0)
    fprintf('The R(Z)=X and I(Z)=Y are uncorrelated\n');
end
c1=corr(X,2*Y);
%c2=cov(X,2*Y);
if (round(c1)==1 || round(c1)==-1)
    fprintf('The R(Z)=X and I(Z)=2*Y are correlated with r = %d\n',c1);
    elseif (round(c1)==0)
    fprintf('The R(Z)=X and I(Z)=2*Y are uncorrelated\n');
end
d1=corr(X,W);
%d2 = cov(X,W);
if (round(d1)==1 || round(d1)==-1)
    fprintf('The R(Z)=X and I(Z)=W are correlated with r = %d\n',d1);
    elseif (round(d1)==0)
    fprintf('The R(Z)=X and I(Z)=W are uncorrelated\n');
end

%%%%%%Checking if sum is Real Gaussian Rv %%%%%%
k1 = X+X; m1 = mean(k1); s1 = std(k1);
k2 = X+Y; m2 = mean(k2); s2 = std(k2);
k3 = X+2*Y; m3 = mean(k3); s3 = std(k3);
k4 = X+W; m4 = mean(k4); s4 = std(k4);

subplot(2,2,1); hist(k1); title('X+X');
subplot(2,2,2);  hist(k2); title('X+Y');
subplot(2,2,3);  hist(k3); title('X+2*Y');
subplot(2,2,4);  hist(k4); title('X+W');

fprintf('\nPaused. Press enter to continue.\n');
fprintf('\n');
pause;

%% Askhsh 4
clc; clear; close all;

n=1000;
k=1000; %z
a1 = rand(n,1);
a2 = rand(n,1);

X1 = -2 + (2+2).*rand(n,1);
m1=mean(X1); %should be 0
var1=var(X1);

X2 = -2 + (2+2).*rand(n,1);
m2 = mean(X2); %should be 0
var2=var(X2);
 
N = -2 + (2+2).*rand(n,n);
I = eye(n,n);
mN = mean(N); %should be 0 nx1
varNs=mean(var(N));	%should be single value
varN=varNs*I;	%true var of N is NxN matrix with varNs on the diagonal

Y=a1*X1'+a2*X2'+N;
mY=mean(Y);
varY=var(Y);

erwthma1=a1*a1'*var1+a2*a2'*var2+varN;

Ry=0;
for i=1:n
    Ry=Ry+Y(:,i)*Y(:,i)';
end
Ry2=Ry/n;

%erotima a, mesh time twn megistwn apoklisewn, prepei na nai mikrh
diafora=mean(max(abs(erwthma1-Ry2)));
fprintf('(erotima a) Difference of theoritical Ry and estimated: %d\n',diafora);																		

%erotima b
xx=rand(n,1);
deksi=mean(var(N))*norm(xx)^2;
aristero=xx'*Ry2*xx;
if (aristero>deksi)
    fprintf('\n(erotima b) isxuei')
end

%erotima c
e=eig(erwthma1);
if (isreal(e) && min(e>0)==1)
    fprintf('\n(erotima c) Eigvalues are real and positive')
end

%erotima d
%megisth timh afaireshs tou varN apo ta n-2 eigvals
fprintf('\n(erotima d)');
d_diff=max((e(1:n-2,:))-varNs*ones(n-2,1))

%z
%Ry(=Z) as Y*Y' and compared to (a)
Z=1/k*Y*Y';
fprintf('\n(erotima z)');
z_diff=max(max(erwthma1-Z))

%h
[U,S,V]=svd(Y);
%h start
h=1/n*U*S*V'*(U*S*V')';
%h end
h_us=U*1/n*S^2*U';
fprintf('\n(erotima h)');
h_diff=max(max(h-h_us))

%theta
fprintf('\n(erotima z)');
detZ=det(Z)
det_temp=1;
for i=i:n
        det_temp=det_temp*S(i,i);
end
det_deksi=det_temp^2/k


%fprintf('\nProgram paused. Press enter to continue.\n');
%pause;



%% FUNCTIONS
% Input: matrix A(m,n) with rank=n, m<n
%Output: matrix Q(m,rank) unitary so that Q'*Q=I
% and matrix R(rank,n) upper triagonal so that A=Q*R
function [Q,R]=myQR(A)
[m,n]=size(A);
r=rank(A);
%1
R(1,1)=norm(A(:,1));
Q(:,1)=A(:,1)/R(1,1);
%n
for z=2:n
    num = A(:,z);
    for i=1:z-1
        %matrix R, w/o z,z elements
        R(i,z)=Q(:,i)'*A(:,z);
        num = num - ( R(i,z)*Q(:,i) );
    end
    R(z,z)=norm(num);
    if (R(z,z) ~= 0)
        %matrix Q
        Q(:,z) = num / R(z,z);
    end
end
%resize
Q=Q(:,1:r);
R=R(1:r,:);
end