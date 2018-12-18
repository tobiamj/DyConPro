function [bport]=dcp_bpfilt_nuisregs(x,Fs,dur_secs,fbot,ftop)

%{
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018

Input:
1. x is the number of time points in the signal
2. Fs is sampling frequency
3. dur_secs

Not Implemented:
1. bandpass 1=bandpass filter; 2=stopband filter

%}

% Fs=.5;x=225;fbot=.01;ftop=.1;dur_secs=size(x,1)/Fs;

siglen=1:(Fs*dur_secs);
bandpass=1;amp=1;bias=0;
filtnum=1/(x*2);
freqvec=filtnum:filtnum:Fs/2;
skeep=find(freqvec>fbot & freqvec<ftop);
sremove=[find(freqvec<fbot) find(freqvec>ftop)];

for loop1=1:length(freqvec)
    sigsin(:,loop1)=amp.*sin(2*(pi).*siglen.*(freqvec(loop1).*2))+bias;
    sigcos(:,loop1)=amp.*cos(2*(pi).*siglen.*(freqvec(loop1).*2))+bias;
end

if bandpass==1
    S=sigsin(:,sremove); % these regressors will be filtered out as a bandpass filter
    C=sigcos(:,sremove);
elseif bandpass==2
    S=sigsin(:,skeep); % these regressors will be filtered out as a notch/stopband filter
    C=sigcos(:,skeep);    
end

bport=[S C];




