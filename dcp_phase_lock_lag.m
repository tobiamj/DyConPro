function [plv,pli,MI]=dcp_phase_lock_lag(x,y)

% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Computes plv and pli between a single time series and all other time
% series channels
% 
% Inputs:
% 1. x is time x channel matrix (can be single column or row)
% 2. y is single time series vector
% 
% Outputs:
% 1. plv is phase locking value
% 2. pli is phase lag index
% 3. MI is cfc modulation index; assumes x is low_sig and y is high_sig amp.
% 
 
phis_x=angle(hilbert(x));

% Computes static plv and pli
phis_y=angle(hilbert(y));
delta_phi=bsxfun(@minus,phis_y,phis_x);
plv=(abs(sum(cos(delta_phi),1)))/size(x,1);
pli=(abs(sum(sign(delta_phi),1)))/size(x,1);

% CFC modulation index (MI)
ab=zscore(abs(hilbert(y))); % amp of hfo
P=ab.*exp(j.*phis_x); % complex signal
MI=abs(mean(P)); % modulation index

