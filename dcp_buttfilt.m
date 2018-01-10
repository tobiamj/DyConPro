function fd=dcp_buttfilt(x,order,LF,HF,fs)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

%  1. input x is time series data
%  2. order is filter order; default = 6
%  3. LF,HF are band stops
%  4. fs is sampling frequency in Hz; e.g., 2 sec TR = .5 Hz
%  5. output y is filtered time series; x is original data

if isempty(order)==1 || exist('order','var')==0
    order=6;
end
if isempty(fs)==1 || exist('fs','var')==0
    fs=1;
end
if isempty(LF)==1 || exist('LF','var')==0
    LF=.25*(fs/2);
end
if isempty(HF)==1 || exist('HF','var')==0
    HF=.75*(fs/2);
end

[b,a]=butter(order,[LF HF]/(fs/2),'bandpass');
fd=filtfilt(b,a,x);


end