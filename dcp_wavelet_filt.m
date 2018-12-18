function [Wf,freqs]=dcp_wavelet_filt(X,Fs)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018

if ~exist('Fs','var')
    Fs=.5;
end

% cardiac_f=[];
% resp_f=[];

[~,rois]=size(X);

for loop1=1:rois
    [wt,freqs]=cwt(X(:,loop1),Fs);
    Wf(:,:,loop1)=wt;
end


end