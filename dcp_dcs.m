function [dcs,ifc]=dcp_dcs(x)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.0 release 10/18/2017

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
% wb=waitbar(0,'Calculating DCS & iFC');
% steps=col;
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
%     waitbar(loop1/steps,wb);
end
% close(wb)

