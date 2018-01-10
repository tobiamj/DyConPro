function zout=dcp_zscore_thresh(x)%,thr)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% Apply zscore normalization on a tensor, matrix or vector
% zout is zscored input, optionally thresholded (not yet!)
% Output is not zscored to preserve values from a non-normal distribution
%   Instead, the values are normalized by their standard deviation

DIMS=ndims(x);

if DIMS==1
    zout=zscore(x);
end
if DIMS==2
    sizediag=length(diag(x));
%     if ~isequal(norm(x),1)
        x=dcp_2Dto1D(x);
        x=zscore(x);
%     end
%     if isequal(round(norm(x),15),1)
% % %         x=dcp_2Dto1D(x);
% % %         x=x./std(x);
%     end
    zout=dcp_1Dto2D(x);%,sizediag);
end
if DIMS==3
    sizediag=length(diag(squeeze(x(1,:,:))));
    [tt,~,~]=size(x);
    for loop1=1:tt
        xx=dcp_2Dto1D(squeeze(x(loop1,:,:))); 
        xx=zscore(xx);
        zout(loop1,:,:)=dcp_1Dto2D(xx,sizediag);
    end
end

% if ~isempty(thr)==1
%     zout(zout<thr)=0;
% end

end