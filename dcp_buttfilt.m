function fd=dcp_buttfilt(x,order,LF,HF,Fs,ftype)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018

% NOW WITH NOTCH FILTERING !! added 12/12/2018 MJT
% 
%  1. input x is time series data
%  2. order is filter order; default = 6
%  3. LF,HF are band stops
%  4. fs is sampling frequency in Hz; e.g., 2 sec TR = .5 Hz
%  5. ftype is a string; 'bandpass' or 'stop'
%  6. output y is filtered time series; x is original data

% clear,clc
% x=randn(1,575);
% order=3;
% LF=.1;
% HF=.14;
% fs=.5;
% ftype='stop';

if isempty(order) || ~exist('order','var')
    order=4;
end
if isempty(Fs) || ~exist('Fs','var')
    Fs=1;
end
if isempty(LF) || ~exist('LF','var')
    LF=.25*(Fs/2);
end
if isempty(HF) || ~exist('HF','var')
    HF=.75*(Fs/2);
end
if ~exist('ftype','var')
    ftype='bandpass';
end

[b,a]=butter(order,[LF HF]./(Fs/2),ftype);
fd=filtfilt(b,a,x);

end

