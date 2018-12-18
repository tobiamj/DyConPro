function [f_alias]=dcp_freq_alias(Fr,Fs)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Inputs:
% 1. Fr=.3; % frequency of (respiration) signal
% 2. Fs=.5; % sampling rate of data
%
% Output:
% 1. f_alias % the center frequency of the Fr aliased into a lower sampling rate
%

Ny=Fs/2; % compute nyquist freq

% f_alias=Ny-mod(Fr,Ny);
f_alias=abs(mod(Fr+Ny,Fs)-Ny);

end
