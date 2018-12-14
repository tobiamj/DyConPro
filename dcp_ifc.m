function [ifc]=dcp_ifc(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018

[row,col]=size(x);
[x1,xfr,~]=dcp_mirror_pad(x);
% ampx=abs(hilbert(x1));
% ampx=ampx(xfr+1:xfr+row,:);
phx=angle(hilbert(x1));
phx=phx(xfr+1:xfr+row,:);
ifc=zeros(row,col,col);
% ampmod=zeros(row,col,col);
% ampx=ampx./norm(ampx);
% phimod=zeros(row,col,col);
% phx=phx./norm(phx);
for slice=1:col
    ifc(:,:,slice)=bsxfun(@minus,phx(:,slice),phx);
%     ifca=bsxfun(@minus,phx(:,slice),phx);
%     ifcb=bsxfun(@plus,phx(:,slice),phx);
%     amp1=bsxfun(@minus,ampx(:,slice),ampx);
%     amp2=bsxfun(@plus,ampx(:,slice),ampx);
%     ampmod(:,:,slice)=1-(abs(amp1)./amp2);
%     phimod(:,:,slice)=1-((ifca)./ifcb);
end
ifc=cos(ifc);

end