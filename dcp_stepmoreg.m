function [BRIK_den,outstrc]=dcp_stepmoreg(Opts)

%
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% function [BRIK_den,outstrc]=dcp_stepmoreg(Opts)
% 
% This function requires data, a brain mask and motion params file.
% 
% Inputs:
% Opts.data=; %path to data file
% Opts.mask=; %path to mask file
% Opts.fwhm_xyz=; %path to file of FWHM_xyz time series output from afni's 3dLocalstat -FWHM
% Opts.moparams=; %path to motion params file
% Opts.brikout=1; %write out a BRIK file? 0=no, 1=yes
% Opts.prefix='testsub'; %name of output BRIK file
% Opts.framewised=1; %compute framewise displacement? 0=no, 1=yes
% Opts.morexs=; %additional params to regress out from data; path to textfile name
% Opts.F24=1; %compute Friston24? 0=no, 1=yes
% Opts.fd_thr=.25; %threshold for framewise displacement scrubbing
% Opts.nodenoise=0; %Only output motion metrics? 0=no, 1=yes
% Opts.view='+orig'; %can be +tlrc or +orig
% Opts.gsr=0; %perform GSR 0=no, 1=yes; It will be expanded with derivs. & squares
% Opts.edgemask; %path to edgemask file
% Opts.edgenum=24; % number of edgemask regressors to keep; default=24
% Opts.remean=0; %Re-mean the data after denoising? 0=no, 1=yes
% Opts.xvarsnorm=0; %Unit norm Xvars; 0=no, 1=yes
% Opts.swthr=0; %apply bonferonni to stepwise threshold; 0=no, 1=yes
% Opts.nostep=0; %perform regular regression instead of stepwise; 0=no, 1=yes
% Opts.tissuemask=''; %path to mask file with CSF=1, GM=2, WM=3;
% Opts.voxcorr=0; %compute voxel-level pairwise correlations to get RMS global correlations; 0=no, 1=yes; will add 30 minutes to your processing time
% Opts.compregs=0; % take first N eigenvectors from svd of Xvars matrix; 0=no, N=yes & how many (N is a number > 0)
% Opts.Fs=; % sampling frequency in Hertz .5 Hz is 2 sec TR
% Opts.fbands=; % fbands is a vector of [fbot ftop] to keep/remove;
%   NOTE: addition of bport to Xvars will occur after all other Xvars have been
%   prepared (i.e., bport is not included in compregs, but is concatenated to
%   the end of the compregs Xvars mat)
%
% 
% Not yet implemented in this function:
% Opts.eigvecdegree=1; % compute eigenvector centrality instead of global correlations
% Opts.nomoregs=0; % do NOT use motion regressors; 0=no (use them), 1=yes (don't use them)
% Opts.maskwmval=3;
% Opts.maskcsfval=1;
% Opts.maskgmval=2;
% Opts.filtxvars=0; %bandpass Xvars & data prior to regression? 0=no, 1=yes
% Opts.physiosigs=; %physiological nuisance regressors
% 
% Outputs:
% 1. _denoised fmri data; _xbeta (beta, R2, adjrsq, t, 1-p); _regorder; _fullmodel (F, 1-p, #df used); _noise_model
% 2. outstrc is a structure with other useful stuff
% 
% Notes:
% 1. Xvars order is: gsr, expanded motion params, expanded tissue regressors (optional), edge
%     regressors (assuming all options are ==1), additional regressors (if any), bandpass filter
%     regressors & fwhm_xyz are the last three regressors to be catenated
%     to Xvars
% 2. Subtract orig data - denoised data to obtain the 'fitted noise signals'
%     NEW !!! Noisemodel is now an automatic output.
% 3. Entering a tissuemask will generate 1 expanded WM and 1 expanded CSF regressors (8 additional regressors)
% 4. Opts.fwhm_xyz could actually be a volume of any voxel specific regressors (No, not implemented to do this)
% 


if ~exist('BrikLoad.m','file') || ~exist('WriteBrik.m','file') || ~exist('afni_matlab','dir')
    error('This function requires the afni_matlab toolbox for the BrikLoad and WriteBrik functions.');
end
if isfield(Opts,'nostep') && Opts.nostep==1
    warning('off', 'stats:pvaluedw:ExactUnavailable');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pathstr,filename,extens]=fileparts(Opts.data);
if strcmp(extens,'.gz')
    filename=filename(1:end-4);
end
[~,BRIK,HEAD,~]=BrikLoad(Opts.data);
[~,MASK,~,~]=BrikLoad(Opts.mask);
motion=load(Opts.moparams);
[xd,yd,zd,td]=size(BRIK);
if isfield(Opts,'tissuemask')
    [~,TMASK,~,~]=BrikLoad(Opts.tissuemask); 
%     GM=TMASK;GM(GM~=2)=0;GM(GM>0)=1;
%     WM=TMASK;WM(WM~=3)=0;WM(WM>0)=1;
%     CSF=TMASK;CSF(CSF~=1)=0;CSF(CSF>0)=1;
%     dilated=zeros(size(GM));
%     for dloop1=1:xd
%         for dloop2=1:yd
%             for dloop3=1:zd
%                 if GM(x,y,z)==1
%                     dilated(x-1:x+1,y-1:y+1,z-1:z+1)=1;
%                 end 
%             end
%         end
%     end
%     WMerode=dilated-WM;WMerode(WMerode<0)=0;
%     CSFerode=dilated-CSF;CSFerode(CSFerode<0)=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'fwhm_xyz') && ~isempty(Opts.fwhm_xyz)
    [~,FWHMxyz,~,~]=BrikLoad(Opts.fwhm_xyz);
    FWHMx=FWHMxyz(:,:,:,1:3:end);
    FWHMy=FWHMxyz(:,:,:,2:3:end);
    FWHMz=FWHMxyz(:,:,:,3:3:end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dmo=diff(motion);
[~,cm]=size(motion);
if cm~=6
    motion=motion';
end
if isfield(Opts,'framewised') && Opts.framewised==1
    fd=zeros(1,size(motion,1)-1);
    for tstep=1:size(motion,1)-1
        fd(:,tstep)=norm(dmo(tstep,:)); % afni style FD=sqrt(sum(motion^2))
    end
else
    fd=[];
end
if isfield(Opts,'F24') && Opts.F24==1
    dmof24=[0 0 0 0 0 0;dmo];
    F24=[motion dmof24 motion.^2 dmof24.^2];
else
    F24=[];
end
if isfield(Opts,'morexs') && ~isempty(Opts.morexs) && Opts.morexs~=0
    xvarsplus=load(Opts.morexs);
else
    xvarsplus=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MASK_rs=reshape(MASK,xd*yd*zd,1);
BRIK_rs=reshape(BRIK,xd*yd*zd,td)';
remeaner=mean(BRIK_rs,1);
this1=find(MASK_rs>=1);
BRIK_mask=BRIK_rs(:,this1);
denoised=zeros(size(BRIK_mask));
censor=zeros(length(fd),1);
censor(fd>Opts.fd_thr)=1;censor=[0 censor']'; % 1 indicates NOT to include this time point (or to spike regress it)
scrubmo=zeros(td,sum(censor));scrubbies=find(censor==1);
for loop1=1:sum(censor)
    scrubmo(scrubbies(loop1),loop1)=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'F24') && Opts.F24==1
    Xvars=F24;
else
    Xvars=motion;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'tissuemask')
% % % % % % % % % % % % % % % %     FIX MASKS TO USE ERODED TISSUE MASKS
    TMASK_rs=reshape(TMASK,xd*yd*zd,1);
    GMmask=find(TMASK_rs==2);GM_sigs=BRIK_rs(:,GMmask);
    WMmask=find(TMASK_rs==3);WM_sigs=BRIK_rs(:,WMmask);
    CSFmask=find(TMASK_rs==1);CSF_sigs=BRIK_rs(:,CSFmask);
    [u,~,~]=svd(zscore(WM_sigs));
    wmvars=u(:,1);wmvars=[wmvars [zeros(1,1);diff(wmvars)] wmvars.^2 [zeros(1,1);diff(wmvars)].^2];
    [u,~,~]=svd(zscore(CSF_sigs));
    csfvars=u(:,1);csfvars=[csfvars [0;diff(csfvars)] csfvars.^2 [0;diff(csfvars)].^2];
    Xvars=[Xvars wmvars csfvars];
else
    wmvars=[];
    csfvars=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'gsr') && Opts.gsr==1
   gsig=mean(BRIK_mask,2);
   dgsig=[0;diff(gsig)];
   Xvars=[gsig dgsig gsig.^2 dgsig.^2 Xvars];
else
    gsig=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'edgemask') && ~isempty(Opts.edgemask)
   [~,emask,~,~]=BrikLoad(Opts.edgemask);
   emask_rs=reshape(emask,1,numel(emask));
   emask_mask=find(emask_rs>=1);
   emaskdat=BRIK_rs(:,emask_mask);
   [u,~,~]=svd(zscore(emaskdat));
   if ~isfield(Opts,'edgenum') || isempty(Opts.edgenum)
       edgenum=24;
   end
   edgevars=u(:,1:Opts.edgenum);
   Xvars=[Xvars edgevars];
else 
    edgevars=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'fwhm_xyz') && ~isempty(Opts.fwhm_xyz)
    Fx_rs=reshape(FWHMx,xd*yd*zd,td)';
    Fy_rs=reshape(FWHMy,xd*yd*zd,td)';
    Fz_rs=reshape(FWHMz,xd*yd*zd,td)';
    Fx_mask=Fx_rs(:,this1);
    Fy_mask=Fy_rs(:,this1);
    Fz_mask=Fz_rs(:,this1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'morexs') && ~isempty(Opts.morexs)
    Xvars=[Xvars xvarsplus];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'swthr') && Opts.swthr==1
    swthr=.05/size(Xvars,2);
else
    swthr=.05;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'nodenoise') && Opts.nodenoise==1
    BRIK_den=[];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'compregs') && ~isempty(Opts.compregs)
    [u,~,~]=svd(zscore(Xvars));
    Xvars=u(:,1:Opts.compregs);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'xvarsnorm') && Opts.xvarsnorm==1    
    for loop1=1:size(Xvars,2)
        Xvars(:,loop1)=Xvars(:,loop1)./norm(Xvars(:,loop1));
    end
end
Xvars=detrend(Xvars); %Xvars & data are detrended
% bandpass filter Xvars & data ?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'Fs') && isfield(Opts,'fbands')
    [bport]=dcp_bpfilt_nuisregs(td,Opts.Fs,td/Opts.Fs,Opts.fbands(1),Opts.fbands(2));
    Xvars=[Xvars bport];
end
% keepxvars=ones(1,size(Xvars,2));
% keepxvars=zeros(1,size(Xvars,2));
% % % % % % % % Xvars=fliplr(Xvars);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stepwise selection of nuisance variables & denoising with linear regression
if isfield(Opts,'nodenoise') && Opts.nodenoise==0
    if ~isfield(Opts,'nostep') || Opts.nostep==0
        if isfield(Opts,'fwhm_xyz') && ~isempty(Opts.fwhm_xyz)
            Xvars=[Xvars Fx_mask(:,1) Fy_mask(:,1) Fz_mask(:,1)];
            betasout=zeros((4*size(Xvars,2))+4,length(this1));
        else
            betasout=zeros((4*size(Xvars,2))+4,length(this1));
        end
        for vox=1:length(this1)
            if isfield(Opts,'fwhm_xyz') && ~isempty(Opts.fwhm_xyz)
               Xvars(:,end-2:end)=[Fx_mask(:,vox) Fy_mask(:,vox) Fz_mask(:,vox)]; 
            end
            tempsig=detrend(BRIK_mask(:,vox));
            [modelb,~,ppp,~,stats,~,history]=stepwisefit(Xvars,tempsig,'display','off','scale','on','penter',swthr); %,'keep',keepxvars);
            denoised(:,vox)=stats.yr; % this is denoised signal
            noiseout(:,vox)=tempsig-stats.yr; % this is the noise model signal
            if ~isempty(history.in)
                fullrsq=1-(stats.rmse^2)/var(tempsig); 
                Adjrsq=fullrsq.*((length(tempsig)-1-stats.df0)/(length(tempsig)-1));
            elseif isempty(history.in)
                fullrsq=0;
                Adjrsq=0;
            end
            modelout(1,vox)=stats.fstat; % this is full model F stat
            modelout(2,vox)=fullrsq;
            modelout(3,vox)=Adjrsq;
            modelout(4,vox)=1-stats.pval; % this is 1-pval for the afni slider
            modelout(5,vox)=stats.df0; % this is the number of df used by the model
            if ~isempty(history.in)
                inlist=find(history.in(end,:)~=0);
                outlist=find(history.in(end,:)==0);
                history.in=history.in(:,inlist);
                rin=sortrows([sum(history.in)' (1:length(sum(history.in)))'],1);
                regin=zeros(1,size(Xvars,2));
                regin(inlist)=rin(:,2); % this is order by which regressors are used by the model
                regorder(:,vox)=regin; % this is order by which regressors are used by the model
                A=zeros(size(modelb));
                adjrsq=abs(diff(history.rmse(rin(:,1)).^2./var(tempsig))); % this is adj. r-squared for each regressor
                adjrsq=[1-(history.rmse(rin(1,1))^2)/var(tempsig) adjrsq];
                for loop1=1:length(adjrsq)
                    A(regorder(:,vox)==loop1)=adjrsq(loop1);
                end
            else
                A=zeros(size(modelb));
                regorder(:,vox)=zeros(size(modelb));
            end
            betasout(1:4:4*length(modelb),vox)=modelb; % this is all betas for the model
            betasout(2:4:4*length(modelb),vox)=A;
            betasout(3:4:4*length(modelb),vox)=stats.TSTAT;
            betasout(4:4:4*length(modelb),vox)=1-ppp; % this is all 1-pvals for all betas in the model
            [length(this1) length(this1)-vox]
        end 
        denoised_full=zeros(td,xd*yd*zd);
        denoised_full(:,this1)=denoised;
        noise_full=zeros(td,xd*yd*zd);
        noise_full(:,this1)=noiseout;
        if isfield(Opts,'remean') && Opts.remean==1
            for loop1=1:size(denoised_full,1)
                denoised_full(loop1,:)=bsxfun(@plus,remeaner,denoised_full(loop1,:));
            end
        end
        BRIK_den=reshape(denoised_full',xd,yd,zd,td);
        betasout(isnan(betasout))=0;betasout(isinf(betasout))=0;
        bfull=zeros(size(betasout,1),xd*yd*zd);bfull(:,this1)=betasout;bfull=reshape(bfull',xd,yd,zd,size(betasout,1));
        modelout(isnan(modelout))=0;modelout(isinf(modelout))=0;
        mfull=zeros(size(modelout,1),xd*yd*zd);mfull(:,this1)=modelout;mfull=reshape(mfull',xd,yd,zd,size(modelout,1));
        regorder(isnan(regorder))=0;regorder(isinf(regorder))=0;
        rfull=zeros(size(regorder,1),xd*yd*zd);rfull(:,this1)=regorder;rfull=reshape(rfull',xd,yd,zd,size(regorder,1)); 
        noise_full=reshape(noise_full',xd,yd,zd,td);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regular linear regression
if isfield(Opts,'nodenoise') && Opts.nodenoise==0
    if isfield(Opts,'nostep') && Opts.nostep==1
        if isfield(Opts,'fwhm_xyz') && ~isempty(Opts.fwhm_xyz)
            Xvars=[Xvars Fx_mask(:,1) Fy_mask(:,1) Fz_mask(:,1)];
            betasout=zeros((4*size(Xvars,2))+4,length(this1));
        else
            betasout=zeros((4*size(Xvars,2))+4,length(this1));
        end
        for vox=1:length(this1)
            if isfield(Opts,'fwhm_xyz') && ~isempty(Opts.fwhm_xyz)
               Xvars(:,end-2:end)=[Fx_mask(:,vox) Fy_mask(:,vox) Fz_mask(:,vox)]; 
            end
            tempsig=detrend(BRIK_mask(:,vox));
            stats=regstats(zscore(tempsig),zscore(Xvars),'linear');
            denoised(:,vox)=stats.r;
            noiseout(:,vox)=stats.yhat;
            modelout(:,vox)=[stats.fstat.f stats.rsquare stats.adjrsquare 1-stats.fstat.pval]'; %F, R2, adj-R2, 1-p
            betasout(1:3:(3*size(Xvars,2))+3,vox)=stats.beta;
            betasout(2:3:(3*size(Xvars,2))+3,vox)=stats.tstat.t;
            betasout(3:3:(3*size(Xvars,2))+3,vox)=1-stats.tstat.pval;
            [length(this1) length(this1)-vox]
        end
        denoised_full=zeros(td,xd*yd*zd);
        denoised_full(:,this1)=denoised;
        noise_full=zeros(td,xd*yd*zd);
        noise_full(:,this1)=noiseout;
        if isfield(Opts,'remean') && Opts.remean==1
            for loop1=1:size(denoised_full,1)
                denoised_full(loop1,:)=bsxfun(@plus,remeaner,denoised_full(loop1,:));
            end
        end
        BRIK_den=reshape(denoised_full',xd,yd,zd,td);  
        betasout(isnan(betasout))=0;betasout(isinf(betasout))=0;
        bfull=zeros(size(betasout,1),xd*yd*zd);bfull(:,this1)=betasout;
        bfull=reshape(bfull',xd,yd,zd,size(betasout,1));
        modelout(isnan(modelout))=0;modelout(isinf(modelout))=0;
        mfull=zeros(size(modelout,1),xd*yd*zd);mfull(:,this1)=modelout;
        mfull=reshape(mfull',xd,yd,zd,size(modelout,1));
        noise_full=reshape(noise_full',xd,yd,zd,td);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE CARPET PLOTS: Pre-/post denoising, fd, noisemodel
if isfield(Opts,'nodenoise') && Opts.nodenoise==0
    cplots=figure('visible','off');
    subplot(4,1,1);imagesc(zscore(BRIK_mask)');colormap 'gray';caxis([-5 5]);colorbar('westoutside')
    subplot(4,1,2);imagesc(zscore(denoised)');colormap 'gray';caxis([-5 5]);colorbar('westoutside')
    subplot(4,1,3);plot(ones(length([0 fd])).*Opts.fd_thr,'k--');hold on;plot([0 fd],'k','LineWidth',2);colorbar('westoutside')
    subplot(4,1,4);imagesc(zscore(noiseout)');colormap 'gray';caxis([-5 5]);colorbar('westoutside')
    saveas(cplots,[Opts.prefix,'_carpet_plots.jpg'],'jpg')
    clear cplots
    close all
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DVARS=[];DVARSQC=[];DVARSMC=[];
if isfield(Opts,'nodenoise') && Opts.nodenoise==0
    DVARS=diff(rms(denoised,2));
    [rr,pp]=corr(fd',DVARS);
    DVARSQC=[rr pp];
    [rr,pp]=corr(fd',diff(rms(BRIK_mask,2)));
    DVARSMC=[rr pp];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmu_denoised=[];rmu_noised=[];
meanr={};nicer={};
mean_denoised=[];mean_noised=[];
if isfield(Opts,'nodenoise') && Opts.nodenoise==0 && isfield(Opts,'voxcorr') && Opts.voxcorr==1
    for loop1=1:size(denoised,2)-1  
        clist=loop1+1:1:size(denoised,2);
        [rr,~]=corr(denoised(:,loop1),denoised(:,clist));
        meanr{loop1}=rr;
        [rr,~]=corr(BRIK_mask(:,loop1),BRIK_mask(:,clist));
        nicer{loop1}=rr;
        length(this1)-loop1
    end
% Fast RMS of voxel-level pairwise correlations
    sizex=(((size(denoised,2)-1)^2)-(size(denoised,2)-1))/2;
    THIS=[];THAT=[];
    for loop1=1:size(denoised,2)-1
       THIS(loop1)=sum(abs(meanr{loop1}).^2);
       THAT(loop1)=sum(meanr{loop1});
    end
    rmu_denoised=sqrt(sum(THIS)/sizex);
    mean_denoised=sum(THAT)/sizex;
    THIS=[];THAT=[];
    for loop1=1:size(denoised,2)-1
       THIS(loop1)=sum(abs(nicer{loop1}).^2);
       THAT(loop1)=sum(nicer{loop1});
    end
    rmu_noised=sqrt(sum(THIS)/sizex);
    mean_noised=sum(THAT)/sizex;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outstrc.inputs=Opts;
outstrc.fd=fd;
outstrc.fd_avg=mean(fd);
outstrc.globalsig=gsig;
outstrc.censor=censor; % 1 means censor it (unlike afni censor files)
outstrc.moscrub=scrubmo; % spike regressors for censored TRs
outstrc.F24=F24; % Friston24 motion params expansion
outstrc.edgevars=edgevars;
outstrc.wmvars=wmvars;
outstrc.csfvars=csfvars;
outstrc.xvars=Xvars; % all regressors that were entered into the model (not necessarily retained)
outstrc.DVARS=[diff(rms(BRIK_mask,2)) DVARS];
outstrc.DVARSQC=[DVARSMC;DVARSQC];
outstrc.rms_global_corr=[rmu_noised rmu_denoised];
outstrc.avg_global_corr=[mean_noised mean_denoised];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(Opts,'brikout') && Opts.brikout==1 && Opts.nodenoise==0
    if isfield(Opts,'nostep') && Opts.nostep==1
        suffix_cell={'_denoised';'_fullmodel';'_xbeta';'_noisemodel'};
        data_cell={'BRIK_den';'mfull';'bfull';'noise_full'};
        BRICK_LABS=cell(4,1);
        for loop2=1:4
            if loop2==1
                for loop1=1:td
                    BRICK_LABS{loop2}=[BRICK_LABS{loop2} ['#',int2str(loop1),'~']];
                end
            end
            if loop2==2
                for loop1=1:size(betasout,1)/4
                    BRICK_LABS{loop2}=[BRICK_LABS{loop2} ['#Full_F_',int2str(loop1),'~#Rsq_',int2str(loop1),'~#AdjRsq_',int2str(loop1),'~#1-pval_',int2str(loop1),'~']];
                end
            end
            if loop2==3
                for loop1=1:size(betasout,1)/3
                    BRICK_LABS{loop2}=[BRICK_LABS{loop2} ['#beta_',int2str(loop1),'~#Rsq_',int2str(loop1),'~#1-pval_',int2str(loop1),'~']];
                end
            end
            if loop2==4
                for loop1=1:td
                    BRICK_LABS{loop2}=[BRICK_LABS{loop2} ['#',int2str(loop1),'~']];
                end
            end
        end
    else
        suffix_cell={'_denoised';'_fullmodel';'_xbeta';'_regorder';'_noise_model'};
        data_cell={'BRIK_den';'mfull';'bfull';'rfull';'noise_full'};
        BRICK_LABS=cell(5,1);
        for loop2=1:5
            if loop2==1
                for loop1=1:td
                    BRICK_LABS{loop2}=[BRICK_LABS{loop2} ['#',int2str(loop1),'~']];
                end
            end
            if loop2==2
                for loop1=1:size(betasout,1)/5
                    BRICK_LABS{loop2}=[BRICK_LABS{loop2} ['#Full_F_',int2str(loop1),'~#FullRsq_',int2str(loop1),'~#AdjRsq_',int2str(loop1),'~#1-pval_',int2str(loop1),'~#inmodel_',int2str(loop1),'~']];
                end
            end
            if loop2==3
                for loop1=1:size(betasout,1)/4
                    BRICK_LABS{loop2}=[BRICK_LABS{loop2} ['#beta_',int2str(loop1),'~#Rsq_',int2str(loop1),'~#tstat_',int2str(loop1),'~#1-pval_',int2str(loop1),'~']];
                end
            end
            if loop2==4
                for loop1=1:size(betasout,1)
                    BRICK_LABS{loop2}=[BRICK_LABS{loop2} ['#regorder_',int2str(loop1),'~']];
                end
            end
            if loop2==5
                for loop1=1:td
                    BRICK_LABS{loop2}=[BRICK_LABS{loop2} ['#',int2str(loop1),'~']];
                end
            end
        end
    end
    for loop1=1:length(data_cell)
        HEAD.BRICK_LABS=BRICK_LABS{loop1};
        Optden.Scale=3;
        Optden.Prefix=[Opts.prefix,suffix_cell{loop1}];
        Optden.View=Opts.view;
        Optden.NoCheck=1;
        Optden.verbose=[];
        Optden.AppendHistory=[];
        Optden.Slices=[];
        Optden.Frames=[];
        Optden.Overwrite='y';
        Optden.AdjustHeader='y';
        [~,~,~]=WriteBrik(eval(data_cell{loop1}),HEAD,Optden);
    end    
end
% Need to include code for saving as mat files if not BRIK/HEAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
