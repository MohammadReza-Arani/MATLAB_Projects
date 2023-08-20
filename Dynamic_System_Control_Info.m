%% azmayesh 4
%Initial Information around Dynamic Systems
clear; clc;
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
num=input("input num[] as a vector  for e.g. [1 0 2 1 ] for s^3+2*s+1  :");%soorate kasr
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
denum=input("input denum[] as a vector ");% Denum

H=tf(num,denum)% Show the System Transfer Function
s=tf('s');

%%

[z,p2,k2]=tf2zp(num,denum)% Finding Zeros and Poles
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
disp("Poles are represenetd as p2 matrix above : ")
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
disp("z vector is the zeros of the system ")
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
disp("K2 is the system gain")
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")


%%
% For Step Response we multiply our system by the term : 1/s whcih is the
% laplace transform of the step signal
%to to that we need to add a 0 at the end of the denum term which shifts
%all the powers of s in the denum by 1
denum2=[denum,0];
disp("Step response of our  system is :")
H2=tf(num,denum2)
%%

% num=[1 1];
% denum2=[1 2+1 2];
[r,p,k]=residue(num,denum2)% Finding System expansion for step response:
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
disp("Poles are represenetd as p matrix above : ")
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
disp("r vector is the residuals for each pole")
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
disp("K is the polynomial term of num/denum after expansion")
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

disp("The fraction expansion contains :")
disp(length(p)+length(k))
disp("Terms as :")
for i=1:length(p)
    disp("           *****************")
    disp("Term "+ i +"  :")
    Hexpanded=tf(r(i),[1,-p(i)])
    clear Hexpanded
end
if(length(k)~=0)
    disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
disp("The last Term (K(s)) is : ")
K=tf(k,[1])
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
end
%% Part 2
G1=tf([-1],[1 0]);
G2=tf(1,[1 3]);
G3=tf(2,[1 9 8]);
G4=tf(2,1);
G5=tf(1,[1 0]);
G6=tf(4,[1 2]);

GT1=feedback(G3,G2);%2 feedbacks in the forward path
GT2=feedback(G5,G6);
GT4=series(GT2,GT1);

GT3=series(G4,G6);% in the feed back


Gforward=series(G2,GT4);% forward path gain
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
disp("The Total transfer function is: ")

Htotal=feedback(Gforward,GT3)% total system with last feedback
disp("           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")

