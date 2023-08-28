clc;
clear;
[x,FS]=audioread('AliLikesToWrite.wav');
sign_x=sign(x);
% sound(sign_x,FS);
Fourierfunc(1,sign_x,FS,'sign_x');
Fourierfunc(1,x,FS,'original');


