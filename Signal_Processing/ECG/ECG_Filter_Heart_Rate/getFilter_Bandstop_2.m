function Hd = getFilter_Bandstop_2
%GETFILTER_BANDSTOP_2 Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.12 and DSP System Toolbox 9.14.
% Generated on: 09-Jan-2023 19:16:32

% Butterworth Bandstop filter designed using FDESIGN.BANDSTOP.

% All frequency values are in Hz.
Fs = 500;  % Sampling Frequency

N   = 70;  % Order
Fc1 = 45;  % First Cutoff Frequency
Fc2 = 55;  % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandstop('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');

% [EOF]
end