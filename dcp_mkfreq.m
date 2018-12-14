function [z,pts]=dcp_mkfreq(ts,fbot,ftop,sr)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% create a vector of frequencies from fbot to ftop Hz
% output vector z length is specified by ts (length of input time series)
% sr is sampling rate in Hz (i.e., same as Fs in other programs)

if length(ts)==1
    z=(fbot:((ftop-fbot)/(ts-1)):ftop);
else
    z=(fbot:((ftop-fbot)/(length(ts)-1)):ftop);
end


if exist('sr','var')
    ts=ts-mean(ts);
    [pts,~]=periodogram(ts,[],z,sr);
%     if exist('plotyes','var')==1 && plotyes == 1
%         figure();plot(w,pts) %PSD
%         [pxp,wp]=periodogram(ts,[],z,sr,'power');
%         figure(20);plot(w,pts/sum(pts)) %POWER i.e., normalized (variance sums to 1)
%     end
end

end