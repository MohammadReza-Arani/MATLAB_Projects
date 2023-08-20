function Hd = getFilter_BandPass_2
%GETFILTER_BANDPASS_2 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.12 and Signal Processing Toolbox 9.0.
% Generated on: 29-Jan-2023 22:55:30

% Elliptic Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 100;  % Sampling Frequency

N      = 50;  % Order
Fpass1 = 0.5;   % First Passband Frequency
Fpass2 = 45;  % Second Passband Frequency
Apass  = 1;   % Passband Ripple (dB)
Astop  = 80;  % Stopband Attenuation (dB)

% Construct an FDESIGN object and call its ELLIP method.
h  = fdesign.bandpass('N,Fp1,Fp2,Ast1,Ap,Ast2', N, Fpass1, Fpass2, ...
                      Astop, Apass, Astop, Fs);
Hd = design(h, 'ellip');

% [EOF]
