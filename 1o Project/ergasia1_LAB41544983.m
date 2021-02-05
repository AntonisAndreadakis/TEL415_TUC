%TEL415 - ERGASIA 1
%ARI8MOS OMADAS: LAB41544983
%Matsatsos Ioannis 2013030148
%Andreadakis Antonis 2013030059
clc; clear all; close all;

%% Askhsh 1
A=[1 0 0 3 3 3; 1 0 0 3 3 3; 1 2 2 0 0 0];
B=[-1 -1; 0 0; 0 0; -1 -2; -1 -2; -1 -2];
C=A*B

A1= [1 0 0; 1 0 0]; A2=3*ones(2,3); A3=[1 2 2]; A4=zeros(1,3);
B1=[-1 -1; 0 0; 0 0]; B2=[-1 -2; -1 -2; -1 -2];
BlockA=[A1 A2; A3 A4]; BlockB=[B1; B2];
BlockC=BlockA*BlockB

disp('^Askhsh 1 - Press a key to run next excercise');
pause
clc; clear all;

%% Askhsh 2
% a)
A_2a=[1 1; 1 2]
A_2a_inv=inv(A_2a)
I=A_2a_inv*A_2a      %yes

% b)
C= [1 1; 1 2]
B=[C C; C 2*C]
B_inv=inv(B)
is_I=B*B_inv      %yes

disp('^Askhsh 2 - Press a key to run next excercise');
pause
clc; clear all;

%% Askhsh 3
A=[1 1; 1 2];
Ai=inv(A);
B=[1 3; 1 4];
Bi=inv(B);
Si=inv(A+B);
inv(Ai+Bi)
A*Si*B
B*Si*A
% all 3 equal

disp('^Askhsh 3 - Press a key to run next excercise');
pause
clc; clear all;

%% Askhsh 4
A=[2 3 1; 1 2 1; 3 3 3];
B=[2 3 1; 4 2 1; 3 3 3];
C=zeros(3,3);
C(2,1)=3;
c=[0; 1; 0;];
d=[3; 0; 0;];
D=A+C;
invB=inv(B);
regular_inverse=inv(D)
rank1_update=inv(A)- ( inv(A)*c*d'*inv(A) )/(1+d'*inv(A)*c)
% both equal

disp('^Askhsh 4 - Press a key to run next excercise');
pause
clc; clear all;

%% Askhsh 6
u1=[3; 5];
u2=[1; 1];
v1=[2; 2; 1];
v2=[3; 1; 3];
A=u1*v1';
B=u2*v2';
if (rank(A+B)==rank(A)+rank(B))
    disp('rank(A+B)==rank(A)+rank(B)');
end

disp('^Askhsh 6 - Press a key to run next excercise');
pause
clc; clear all;

%% Askhsh 7
A = [1 1 0; -1 1 2; 1 1 0; -1 1 2]
       
rankA=rank(A)           % r equals dim R(A)
           
basisRA=orth(A)       % num of cols equals r
basisNA=null(A)

disp('---------------erwthma a-----------------')
b=[0; 0; 0; 0]
rref([basisRA b])       % if b exists in R(A)? b does! So Ax=b can be solved!
x0=[1;-1;1]                % one manual initial solution
x=x0+1*basisNA      % generate all solutions
A*x                            %should be b
%confirm that above x solve A'*A*x=A'*b
if(int64(A'*A*x)==int64(A'*b))
    disp('They do! (a)');
end

disp('---------------erwthma b-----------------')
b=[2; 0; 2; 0]
rref([basisRA b])       % if b exists in R(A)? b does! So Ax=b can be solved!
x0=[1;1;0]
x=x0+2*basisNA
A*x                             %should be b
%confirm that above x solve A'*A*x=A'*b
if(int64(A'*A*x)==int64(A'*b))
    disp('They do! (b)');
end

disp('---------------erwthma c-----------------')
b=[0; 2; 2; 0]
rref([basisRA b])       % if b exists in R(A)? b DOESN'T. Only A'*A*x=A'*b can be solved!
% ???
disp('^Askhsh 7');