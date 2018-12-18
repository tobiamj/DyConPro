function [IFC,IFC_std,ifc_dfcmap,ifc_map,outstrc]=dcp_fast_s2b_ifc(seeds,X,mask,bandpass,prefix,outview,suffix,surrogate)

% 
% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% This function will compute the DCS between each voxel from a seed region, 
% and all other voxels in the brain.  It has two outputs: connectivity,
% and connectivity variability

% Usage:
% 1. data is a string with full path and filename of the .BRIK or .nii file to load
% 2. seeds is the voxel location to use as seed region [x y z], or mask index #,
%       or array of seed time courses
% 3. prefix is a string with the name of the output file (no it's not; yes it is)
% 4. bandpass is a [1 2 3] with fbot btop Fs
% 5. mask is a regional parcellation mask with integers identifying index #
% 6. DCS (or IFC) is the output
% 7. DCS_std (or IFC_std) is the TR-wise variability of the seed voxels' DCS
% 8. surrogate is for computing phase randomized confidence intervals; 0=no, 1=yes
% 
% Notes:
% 1. bandpass filter parameters are hard coded
% 2. this function will apply an internally generated mask if a mask is not provided
% 3. Data and seed time courses are zscored column-wise prior to any further calculations
% 

% Load data if it exists
if ischar(X)
    [fd1,~,ext]=fileparts(X);
    cd(fd1)
    if strcmp(ext,'.BRIK')
        if ~exist('BrikLoad.m','file') || ~exist('WriteBrik.m','file')
            error('You need the afni_matlab toolbox to load your data.');
        end
        [~,V,BHEAD,~]=BrikLoad(X);
        [xd,yd,zd,td]=size(V);
        X=reshape(X,xd*yd*zd,td)';
        X=zscore(X);
    elseif strcmp(ext,'.nii')
        if ~exist('spm_read_vols.m','file') || ~exist('spm_vol.m','file')
            error('You gotta get some spm functions to load your data.');
        end
        V=spm_vol(X);
        [X,~]=spm_read_vols(V);
        [xd,yd,zd,td]=size(X);
        X=reshape(X,xd*yd*zd,td)';
        X=zscore(X);
    end
elseif ndims(X)==4
    [xd,yd,zd,td]=size(X);
    X=reshape(X,xd*yd*zd,td)';
    X=zscore(X);
else  
    error('What is wrong with your data?')
end

% Load seedmask
if exist('seeds','var') && ischar(seeds)
    [fs1,~,fext]=fileparts(mask);
    cd(fs1)
    if strcmp(fext,'.BRIK')
        [~,SMASK,BHEAD,~]=BrikLoad(seeds);
        SMASK=reshape(SMASK,1,xd*yd*zd);
    elseif strcmp(fext,'.nii')
        V=spm_vol(seeds);
        [SMASK,~]=spm_read_vols(V);
        SMASK=reshape(SMASK,1,xd*yd*zd);
    end
    maskchar=1;
elseif ndims(seeds)==4 || ndims(seeds)==3
    SMASK=reshape(seeds,1,xd*yd*zd);
    maskchar=1;
else
    SMASK=X(:,:,:,1);
    SMASK(SMASK~=0)=1;
    SMASK=reshape(SMASK,1,xd*yd*zd);
end
cd(fd1)

% Load mask if it exists
if exist('mask','var') && ischar(mask)
    [fm1,~,ext]=fileparts(mask);
    cd(fm1)
    if strcmp(ext,'.BRIK')
        [~,MASK,~,~]=BrikLoad(mask);
        MASK=reshape(MASK,1,xd*yd*zd);
    elseif strcmp(ext,'.nii')
        V=spm_vol(mask);
        [MASK,~]=spm_read_vols(V);
        MASK=reshape(MASK,1,xd*yd*zd);
    end
    maskchar=1;
elseif ndims(mask)==4 || ndims(mask)==3
    MASK=reshape(mask,1,xd*yd*zd);
    maskchar=1;
else
    MASK=X(:,:,:,1);
    MASK(MASK~=0)=1;
    MASK=reshape(MASK,1,xd*yd*zd);
end
cd(fd1)
seedtot=unique(SMASK(SMASK~=0));

% % % % ADAPT FOR MULTIPLE UNIQUE SEEDS; Oh nevermind... One seed region is good
for loop1=1:length(seedtot)
    seedts(:,:,loop1)=X(:,SMASK==seedtot(loop1));
end

% Find info about and apply the mask
% this0=find(MASK==0);
this1=find(MASK~=0);
X=X(:,this1);
voxtot=length(this1);

% bandpass filter the data if desired
if exist('bandpass','var') && length(bandpass)==3
    X=dcp_buttfilt(X,3,bandpass(1),bandpass(2),bandpass(3));
    seedts=dcp_buttfilt(seedts,3,bandpass(1),bandpass(2),bandpass(3));
end

% create phase vectors
AHBRIK=angle(hilbert(X));
% AbsBRIK=abs(hilbert(X));
seeds1=angle(hilbert(seedts));

randal=randi(size(seeds1,1),size(seeds1,1),size(seeds1,2));
randseeds1=seeds1(randal);

% pre-allocate some matrices
cosync=zeros(td,voxtot);
cosync_std=zeros(td,voxtot);

% main ifc/dcs loop
for loop1=1:voxtot  
    sync=cos(bsxfun(@minus,AHBRIK(:,loop1),seeds1));    
    cosync(:,loop1)=mean(sync,2);
    cosync_std(:,loop1)=std(sync,[],2);
    [voxtot voxtot-loop1]
    if exist('surrogate','var') && surrogate==1
        randsync=cos(bsxfun(@minus,AHBRIK(:,loop1),randseeds1));    
        randcosync(:,loop1)=mean(randsync,2);
        randcosync_std(:,loop1)=std(randsync,[],2);
    end
end
% reshape IFCs to match BRIK size
IFC=zeros(td,xd*yd*zd);IFC_std=zeros(td,xd*yd*zd);
IFC(:,this1)=cosync;IFC_std(:,this1)=cosync_std;
IFC=reshape(IFC',xd,yd,zd,td);
IFC_std=reshape(IFC_std',xd,yd,zd,td);
ifc_dfcmap=std(IFC,[],4);
ifc_map=mean(IFC,4);

% ifc_map=ifc_map./max(max(max(ifc_map)));

% determine subject-specific confidence intervals
RNK=sortrows(reshape(randcosync,numel(randcosync),1),1);
mRNK=mean(RNK);
stdRNK=std(RNK,[],1);
zcrit=[1.96 2.326 2.576 3.090 3.719]; % critical Z-scores for 95%, 99%, 99.5%, 99.9% & 99.99% confidence intervals

outstrc.CI95=[mRNK-(zcrit(1)*stdRNK) mRNK+(zcrit(1)*stdRNK) .05];
outstrc.CI99=[mRNK-(zcrit(2)*stdRNK) mRNK+(zcrit(2)*stdRNK) .01];
outstrc.CI995=[mRNK-(zcrit(3)*stdRNK) mRNK+(zcrit(3)*stdRNK) .005];
outstrc.CI999=[mRNK-(zcrit(4)*stdRNK) mRNK+(zcrit(4)*stdRNK) .001];
outstrc.CI9999=[mRNK-(zcrit(5)*stdRNK) mRNK+(zcrit(5)*stdRNK) .0001];

% save outputs as BRIK/HEAD files
Optdcs.Scale=1;
Optdcs.Prefix=[prefix,'_IFC_S2B',suffix];
Optdcs.View=outview;
Optdcs.NoCheck=1;
Optdcs.verbose=[];
Optdcs.AppendHistory=[];
Optdcs.Slices=[];
Optdcs.Frames=[];
Optdcs.Overwrite=[];
Optdcs.AdjustHeader='y';

Optvar.Scale=1;
Optvar.Prefix=[prefix,'_IFC_var',suffix];
Optvar.View=outview;
Optvar.NoCheck=1;
Optvar.verbose=[];
Optvar.AppendHistory=[];
Optvar.Slices=[];
Optvar.Frames=[];
Optvar.Overwrite=[];
Optvar.AdjustHeader='y';

Optvarb.Scale=1;
Optvarb.Prefix=[prefix,'_IFC_dfcmap',suffix];
Optvarb.View=outview;
Optvarb.NoCheck=1;
Optvarb.verbose=[];
Optvarb.AppendHistory=[];
Optvarb.Slices=[];
Optvarb.Frames=[];
Optvarb.Overwrite=[];
Optvarb.AdjustHeader='y';

Optvarc.Scale=1;
Optvarc.Prefix=[prefix,'_IFC_map',suffix];
Optvarc.View=outview;
Optvarc.NoCheck=1;
Optvarc.verbose=[];
Optvarc.AppendHistory=[];
Optvarc.Slices=[];
Optvarc.Frames=[];
Optvarc.Overwrite=[];
Optvarc.AdjustHeader='y';

[~,~,~]=WriteBrik(IFC,BHEAD,Optdcs);
[~,~,~]=WriteBrik(IFC_std,BHEAD,Optvar);
[~,~,~]=WriteBrik(ifc_dfcmap,BHEAD,Optvarb);
[~,~,~]=WriteBrik(ifc_map,BHEAD,Optvarc);


% keyboard
% 
% % cd (pathstr)
% if ~exist('./DCS_S2B','dir')
%     mkdir ('DCS_S2B')
% end
% cd (['./DCS_S2B'])
% 
% end
