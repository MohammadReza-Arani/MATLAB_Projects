clc;
clear;
[y,Fs]=audioread('tel.wav');
window=hamming(100);
noverlap=99;
nfft=100;
HX2=spectrogram(y,window,noverlap,nfft,Fs);
figure
imagesc(abs(HX2));
title('Spectrogram output');
