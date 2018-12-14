% function [outfile]=dcp_nuis_motion_f24(dcpfile,motionfile,options)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% INPUT:
% 1. dcpfile is data file of raw (minimally processed) signals
% 2. motion file is file of six rigid body motion parameters
% 3. options has some options in it for making figure, regressing out 
%     additional nuisance variables and saving filenames

% OUTPUT:
% 1. dcp structure with cleaned data
% 2. files with Friston24 and FD

% generate Friston24 motion parameters and framewise displaement (FD)
for loop1=1:16
    X=load(['subj_',int2str(loop1),'.motion']);
    Xd=[0 0 0 0 0 0;diff(X)];
    Xmo=[X X.^2 Xd Xd.^2];
    fd=sum(abs(Xd),2);
%     save(['subj_',int2str(loop1),'.friston24.txt'],'Xmo','-ascii')
    FDt(:,loop1)=fd;
    FD1(loop1)=mean(fd(1:59));
    FD2(loop1)=mean(fd(60:99));
    FD3(loop1)=mean(fd(100:159));
    FD4(loop1)=mean(fd(160:169));
    FD5(loop1)=mean(fd(170:end));
    for loop2=1:225
        NX(loop2,loop1)=norm(Xd(loop2,:));
    end
end

% generate figure for motion params and FD
% for loop2=1:6 
%     figure();plot(squeeze(diff(X(:,loop2,:))))
% end

% regress out Friston24, an doptionally FD
% for loop2=1:16
%     load(['subj_',int2str(loop2),'_dcp.mat']);
%     load(['subj_',int2str(loop2),'.friston24.mat']);
%     wmsig=load(['subj_',int2str(loop2),'.wmsig.1D']);
%     covars=[Xmo wmsig];
%     for roi=1:200
%         statsout=regstats(dcp.data.ts.raw(:,roi),covars,'linear','r');
%         dcp.data.ts.motion24(:,roi)=statsout.r;
%     end
% %     save(['subj_',int2str(loop2),'_dcp.mat'],'dcp','-mat') 
%     [dcs,~]=dcp_dcs(dcp.data.ts.motion24);
%     dcp.dynfc.dcs=dcs;
%     save(['subj_',int2str(loop2),'_dcp.mat'],'dcp','-mat')
% end






