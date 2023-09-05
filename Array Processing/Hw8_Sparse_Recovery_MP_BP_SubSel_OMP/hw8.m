% be name khoda - Solution of Hw8 - part e
load('hw8.mat')
[M,N]=size(D);
N0=3;
% Linear Programming
f=ones(2*N,1);
Aeq=[D -D];
beq=x;
lb=zeros(2*N,1);
yhat = linprog(f,[],[],Aeq,beq,lb,[]);
splus=yhat(1:N);
sminus=yhat(N+1:end);
sBP=splus-sminus;
posBP=find(abs(sBP)>0.01)';
disp('BP:')
[posBP;sBP(posBP)'] 

    