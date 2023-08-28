clc;
clear;
[y,Fs]=audioread('tel.wav');
 for j=1:7700-50
     Xt(:,j)=y(j:j+50);
 end
 HX=(fft(Xt));
 P=fftshift(HX);
 P2=abs(P);
 figure
 imagesc(P2);
 title('Our Func');
 grid on
