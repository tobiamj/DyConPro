function [dcs,ifc]=dcp_dcs(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018
% 
% Inputs:
% 1. x is a RxC matrix of time series with time down the column and 
%   channels across the columns
% 
% Outputs:
% 1. dcs is the wideband filtered instantaneous time-varying correlation 
%   of two signals. This will be a tensor of size CxCxR
% 2. ifc is the wideband (or narrow if pre-filtered) instantaneous 
%   time-varying correlation of two signals. This will be a tensor of size CxCxR
% 
% Default: output tensors have diagonals set to 0 at each time point
% 
% Notes: you might want to use (abs(dcs)) to get measure of coherence
%   independent of correlation/anticorrelation
% 

rng(1);
[row,col]=size(x);
scaler=mean(range(x))*.1;
noise1=randn(15,1).*scaler;
noise1=repmat(noise1,1,col);
noise2=randn(15,1).*scaler;
noise2=repmat(noise2,1,col);
x1=[noise1;x;noise2];
phx=angle(hilbert(x1));
phx=phx(16:end-15,:);
ifc=zeros(row,col,col);
    
for slice=1:col
    ifc(:,:,slice)=bsxfun(@minus,phx(:,slice),phx);
end
% ifc(ifc>pi)=(2*pi)-ifc(ifc>pi);
ifc=cos(ifc);

dcs=ifc;
wb=waitbar(0,'Calculating DCS & iFC');
steps=col;
for loop1=1:col
    for loop2=1:col
        temp=squeeze(dcs(:,loop1,loop2));
        if loop1~=loop2
        ue=dcp_getspline(temp',[]);
        le=-dcp_getspline(-temp',[]);
        plve=(ue+le)/2;
        dcs(:,loop1,loop2)=plve;
        else
        dcs(:,loop1,loop2)=zeros(1,length(temp));
        end
    end
    waitbar(loop1/steps,wb);
end
close(wb)

