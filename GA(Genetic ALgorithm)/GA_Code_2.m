clear; clc; close all;



f = [350*5,330*3,310*4,280*6,500,450,400,100];

intcon = 1:4;
A = [];
b = [];
Aeq = [5,3,4,6,1,1,1,1;
    5*0.05,3*0.04,4*0.05,6*0.03,0.08,0.07,0.06,0.03;
    5*0.03,3*0.03,4*0.04,6*0.04,0.06,0.07,0.08,0.09];
beq = [25;1.25;1.25];

lb = zeros(8,1);
ub = ones(8,1);
ub(5:end) = Inf; % No upper bound on noninteger variables


[x,fval] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);

nvars = length(Aeq) ;
nonlcon =[];
opts = optimoptions('ga','PlotFcn',@gaplotbestf);

ga_fun = @(x) f(1)*x(1)+f(2)*x(2)+f(3)*x(3)+f(4)*x(4);

[x_ga , fval_ga] = ga(ga_fun,nvars,A,b,Aeq,beq,lb,ub,nonlcon,intcon,opts)




