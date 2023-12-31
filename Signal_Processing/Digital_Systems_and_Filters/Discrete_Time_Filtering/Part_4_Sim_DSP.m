function Hd = Part_4_Sim_DSP
%PART_4_SIM_DSP Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.12 and Signal Processing Toolbox 9.0.
% Generated on: 14-May-2023 23:41:13

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 16000;  % Sampling Frequency

N     = 40;    % Order
Fpass = 2000;  % Passband Frequency
Fstop = 2500;  % Stopband Frequency
Wpass = 1;     % Passband Weight
Wstop = 1;     % Stopband Weight
dens  = 20;    % Density Factor

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop], ...
           {dens});
Hd = dfilt.dffir(b);

% [EOF]
end
