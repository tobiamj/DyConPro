function [roi_tcs]=dcp_roi_timecourse_extract(data,mask,method,filt_pb,dtrnd)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.1 release 12/18/2018
% 
% This function will compute ROI time courses with either of 2 methods.  It
% takes fMRI volumes and an MRI mask parcellation in AFNI format as inputs.
% 
% Inputs:
% 1. data is a string with path/filename+orig.BRIK; data should be fully
%     preprocessed before taking their regional averages (ideally)
% 2. mask is string with path/filename+orig.BRIK with parcellation into N ROIs
% 3. method==1 is average signal in an roi; method==2 is first eigenvector
% 4. filt_pb is filter passband [lf hf fs], e.g., [.01 .15 .5]; THIS WILL
%     BANDPASS FILTER THE SIGNALS BEFORE TAKING THEIR REGIONAL AVERAGE
% 
% Outputs:
% 1. roi_tcs_mean is a Time x ROIs matrix ready for dFC processing
% 
% NOTES:
% 1. Method 1 & 2 will produce virtually identical results.  To make 
%     method 1 equal to method 2 just divide method 1 by the norm of the mean signal
%     for example, roi_tcs_mean(:,1)./norm(roi_tcs_mean(:,1)) == U(:,1)
% 
% 2. Global signal regression is hardcoded; default=0 for do NOT do it !!!
%     If you want to perform global signal regression then this option will compute
%     the global signal over all voxels in the inclusion mask and then regress it 
%     from each voxel after detrending the voxel time series
% 
% 

if ~exist('BrikLoad.m','file') || ~exist('WriteBrik.m','file') || ~exist('afni_matlab','dir')
    error('This function requires the afni_matlab toolbox.');
end

gsr=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndims(data)~=4
    [dpath,dfname,ext]=fileparts(data);
else
    ext=[];
end
if ndims(mask)~=3
    [mpath,mfname,mext]=fileparts(mask);
else
    mext=[];
end

if ~exist('data','var') || isempty(data)
    error('Enter some data and options please')
elseif ~isempty(ext) && strcmp(ext,'.nii') 
    V=spm_vol(data);
    [BRIK,~]=spm_read_vols(V);
elseif ~isempty(ext) && strcmp(ext,'.BRIK') || strcmp(ext,'.gz')
    [~,BRIK,~,~]=BrikLoad(data);
elseif isempty(ext) && ndims(data)==4
    BRIK=data;
end

if ~exist('mask','var') || isempty(mask)
%     [~,MASK,~,~]=BrikLoad('/home/mtobia/BN_Atlas_246_2mm_PhysLearn+tlrc.BRIK');
elseif ~isempty(mext) && strcmp(mext,'.nii') || ~isempty(mext) && strcmp(mext,'.gz')
    if strcmp(mext,'.gz')
        mask=mask(1:end-3);
    end
    Vm=spm_vol(mask);
    [MASK,~]=spm_read_vols(Vm);
elseif ~isempty(mext) && strcmp(mext,'.BRIK')
    cd([filesep,mpath])
    [~,MASK,~,~]=BrikLoad(mask);  
elseif isempty(mext) && ndims(mask)==3
    MASK=mask;
end
cd([filesep,dpath])

if ~exist('method','var') || isempty(method)
    method=1;
end
if ~exist('filt_pb','var') || isempty(filt_pb)
    filt_pb=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xd,yd,zd,td]=size(BRIK);
MASK_rs=reshape(MASK,1,xd*yd*zd);
BRIK_rs=reshape(BRIK,xd*yd*zd,td)';
this0=find(MASK_rs==0);
% this1=find(MASK_rs>0); % Not used in this function
BRIK_mask=BRIK_rs;BRIK_mask(:,this0)=[];
MASK_rs(this0)=[];
voxels=size(BRIK_mask,2);
rois=length(unique(MASK_rs));
roi_tcs=zeros(td,rois);

if gsr==1
    globalsig=mean(BRIK_mask,2);
    Xmoreg=zeros(size(BRIK_mask));
    for loop1=1:size(BRIK_mask,2)
        tempsig=detrend(BRIK_mask(:,loop1)); % detrend linear
        [~,~,Xmoreg(:,loop1),~]=regress(tempsig,globalsig); % global sig regression
    end
    BRIK_mask=Xmoreg;
end

if exist('filt_pb','var')==1 && ~isempty(filt_pb)==1
    [mpad,mfront,~]=dcp_mirror_pad(BRIK_mask);
%     Wp=[filt_pb(1) filt_pb(2)]/(filt_pb(3)/2); 
%     Q=.75;Ws=[Wp(1)*Q Wp(2)/Q];
%     Rp=3;Rs=60;
%     [bw_ord,~]=buttord(Wp,Ws,Rp,Rs);
    bw_ord=4;bpfilt=zeros(size(mpad));
    for vox=1:voxels
        bpfilt(:,vox)=dcp_buttfilt(mpad(:,vox),bw_ord,filt_pb(1),filt_pb(2),filt_pb(3));
    end
    BRIK_mask=bpfilt(mfront+1:end,:);
    BRIK_mask=BRIK_mask(1:td,:);
end

if dtrnd==1
    BRIK_mask=detrend(BRIK_mask);
end

if method==1
% METHOD 1: average signal in each ROI
    for roi=1:rois
        roi_tcs(:,roi)=mean(BRIK_mask(:,MASK_rs==roi),2);
        if exist('filt_pb','var')==1 && ~isempty(filt_pb)==1
            roi_tcs(:,roi)=dcp_buttfilt(roi_tcs(:,roi),bw_ord,filt_pb(1),filt_pb(2),filt_pb(3));            
        end
    end
end

if method==2
% METHOD 2: first eigenvector in each ROI
    for roi=1:rois
        [U,~,~]=svd(BRIK_mask(:,MASK_rs==roi));
        roi_tcs(:,roi)=U(:,1);
        if exist('filt_pb','var')==1 && ~isempty(filt_pb)==1
            roi_tcs(:,roi)=dcp_buttfilt(roi_tcs(:,roi),bw_ord,filt_pb(1),filt_pb(2),filt_pb(3));            
        end
    end
end


end
