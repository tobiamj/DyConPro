function [filt_sigs,outstrc]=dcp_dyadic_buttfilt(X,N,Fs)

% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% Inputs:
% 1. X is time x channel matrix
% 2. N is number of frequency bins from 0 to Nyquist
%     leave this empty, [], to use intrinsic dyadic bands from 0 to Nyquist
%         n.b., this essentially approximates an EMD
% 3. Fs is sampling frequency in Hz, e.g. a TR=2 sec is .5 Hz
% 
% Note:
% 1. adjust f_bot and f_top as needed; they are applied to w_centers, so
%     the actual frequency content may exceed either constraint
% 2. adjust fnum as needed; default fnum=8 (8 freq bands)
%     this will not determine the number of freq bands used in the filter;
%     that is determined by f_bot and f_top constraints on w_centers
% 

f_bot=.001;
f_top=Fs/2;
fnum=10;

sigs=size(X,2);
siglength=size(X,1);
Nyquist=Fs/2;
butt_order=3;
    
if exist('N','var') && ~isempty(N)
    if N<7 
        N=15;
    end
    binstep=Nyquist/N;
    bandrow=0:binstep:Nyquist;bandrow([1 end])=[];
    fbands=[bandrow(1:end-1)' bandrow(2:end)'];
    fbands(1,2)=fbands(1,2)-.0001;
    fbands_full=fbands;
end

% Determine dyadic frequency bands
if ~exist('N','var') || isempty(N)
    fbands=zeros(fnum,2);
    for loop2=1:fnum
        fbands(loop2,:)=[Nyquist/(2^loop2) 2.*(Nyquist/(2^loop2))]; 
    end
    fbands_full=fbands;fbands(1,2)=fbands(1,2)-.0001;
    w_centers=mean(fbands,2);
    centers_keep=find(w_centers>=f_bot & w_centers<=f_top);
    fbands=fbands(centers_keep,:);
end

[Xmp,mf,~]=dcp_mirror_pad(X);

freqs=size(fbands,1);
x_filt_mat=zeros(size(Xmp,1),length(fbands),sigs);
for loop1=1:sigs
    sigX=Xmp(:,loop1);
    for loop2=1:freqs
        x_filt=dcp_buttfilt(sigX,butt_order,fbands(loop2,1),fbands(loop2,2),Fs);
        sigX=sigX-x_filt;
        x_filt_mat(:,loop2,loop1)=x_filt;
    end
end

x_filt_mat=x_filt_mat(mf+1:end,:,:);
x_filt_mat=x_filt_mat(1:siglength,:,:);

% determine if filtered signals are residue
residue=zeros(size(x_filt_mat,2),size(x_filt_mat,3));
for loop1=1:size(x_filt_mat,3)
    for loop2=1:size(x_filt_mat,2)
        [z,~]=dcp_zero_xings(x_filt_mat(:,loop2,loop1));
        if z<3
            residue(loop2,loop1)=1;
        else
            residue(loop2,loop1)=0;
        end
    end
end

filt_sigs=x_filt_mat;
if exist('fbands','var')
    outstrc.fbands=fbands;
end
if exist('w_centers','var')
    outstrc.wcenters=w_centers;
end
if exist('fbands_full','var')
    outstrc.fbands_full=fbands_full;
end
outstrc.residue=residue;


end
