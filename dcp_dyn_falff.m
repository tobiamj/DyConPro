function [dyn_out]=dcp_dyn_falff(Opts)

% Inputs:
% 1. X is input time series data arranged with time down the rows, or it
%   can be a 4D nifti or brik/head file format
% 2. Fs is sampling frequency
% 3. keep_band is the freq bins to keep for computation of alff & falff
% 4. mask is a mask in same format as X
% 
% Output:
% 1. dyn_falff is time-varying falff measure (also malff and alff)
% 
% Notes:
% 1. requires AFNI installation and available on the path
% 2. Opts is structure with thefollowing fields:
%     Opts.X='/home/mtobia/Documents/MATLAB/SMOOTHED.nii';
%     Opts.mask='/home/mtobia/Documents/MATLAB/SMOOTHED.nii';
%     Opts.Fs=.5
%     Opts.wl=15;
%     Opts.overlap=14;
%     Opts.keepband=[.01 .08];
%     Opts.brikout=1; 1=yes; 2=nifti;


% X='/home/mtobia/ScratchMoFolder/HIVCB/RS/193/193-restingstate.feat/rsfmri-prep-nobp-smooth-mni.nii.gz';
% mask='/home/mtobia/BN_Atlas_246_2mm.nii.gz';
% Opts.brikout=2;
% 
% Fs=.5;
% wl=15;
% overlap=14;
% keep_band=[.01 .08];

gsr=1;

X=Opts.X;
mask=Opts.mask;
Fs=Opts.Fs;
wl=Opts.wl;
overlap=Opts.overlap;
keep_band=Opts.keepband;
prefix=Opts.prefix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DATA load
if ~exist('X','var') || isempty(X)
    error('Enter some data and options please')
end
if ischar(X)
    [~,filename,ext]=fileparts(X);
else
    ext=[];
end
if ~isempty(ext) && strcmp(ext,'.nii') 
    V=spm_vol(X);
    [BRIK,~]=spm_read_vols(V);
elseif ~isempty(ext) && strcmp(ext,'.BRIK') || strcmp(ext,'.gz')
    [~,BRIK,HEAD,~]=BrikLoad(X);
elseif isempty(ext) && ndims(X)==4
    BRIK=X;
elseif ismatrix(X)
    BRIK=X;
    Opts.brikout=0;
end
if ~ismatrix(BRIK)
    [xd,yd,zd,td]=size(BRIK);
    BRIK_rs=zscore(reshape(BRIK,xd*yd*zd,td)');
end

%MASK load
if exist('mask','var') || ~isempty(mask)
    if ischar(mask)
        [~,~,mext]=fileparts(mask);
    else
        mext=[];
    end
    if ~isempty(mext) && strcmp(mext,'.nii') 
        V=spm_vol(mask);
        [MASK,~]=spm_read_vols(V);
    elseif ~isempty(mext) && strcmp(mext,'.BRIK') || strcmp(mext,'.gz')
        [~,MASK,~,~]=BrikLoad(mask);
    elseif isempty(mext) && ndims(X)==3
        MASK=mask;
    elseif ismatrix(mask)
        MASK=mask;
    end
end
if ~ismatrix(MASK)
    MASK_rs=reshape(MASK,xd*yd*zd,1);
end

if exist('prefix','var') && ~isempty(prefix) && ischar(prefix)
    filename=prefix;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

this1=find(MASK_rs>=1);
BRIK_mask=BRIK_rs(:,this1);
nsigs=size(BRIK_mask,2);
siglength=size(BRIK_mask,1);

if gsr==1
    BRIK_mask=detrend(BRIK_mask);
    glosig=mean(BRIK_mask,2);
    denoised=zeros(size(BRIK_mask));
    for loop1=1:nsigs
       [~,~,denoised(:,loop1),~]=regress(BRIK_mask(:,loop1),glosig); 
    end
    BRIK_mask=denoised;
end

for loop1=1:nsigs
    
    [x,xfr,~]=dcp_mirror_pad(BRIK_mask(:,loop1));
    [~,F,~,P]=spectrogram(x,wl,overlap,[],Fs);
    P=P(:,xfr+1:xfr+siglength);
    for loop2=1:siglength        
        s1(loop2)=sum(P([F>keep_band(1) & F<keep_band(2)],loop2));
        s2(loop2)=sum(P(:,loop2));
    end
    dfalff(:,loop1)=(s1./s2)';  
    dalff(:,loop1)=s1';
    
    loop1/nsigs
    
end

dmalff=dalff./nanmean(nanmean(dalff,1),2);
dfalff=dfalff./nanmean(nanmean(dfalff,1),2);

dfalff_var=std(dfalff,[],1);
dmalff_var=std(dmalff,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dyn_out.dfalff=dfalff;
dyn_out.dmalff=dmalff;
dyn_out.dfalff_var=dfalff_var;
dyn_out.dmalff_var=dmalff_var;

alffout={'dmalff';'dfalff';'dmalff_var';'dfalff_var'};
for loop1=1:4
    if loop1<3
        d_full=zeros(td,xd*yd*zd);
        d_full(:,this1)=eval(alffout{loop1});
        brikout=reshape(d_full',xd,yd,zd,td);
    end
    if loop1>2
        d_full=zeros(1,xd*yd*zd);
        d_full(:,this1)=eval(alffout{loop1});
        brikout=reshape(d_full',xd,yd,zd,1);
    end
    if isfield(Opts,'brikout') && Opts.brikout>0
        Optden.Scale=0;
        Optden.Prefix=[filename,'_gsr_',alffout{loop1}];
        Optden.View='+tlrc';
        Optden.NoCheck=1;
        Optden.verbose=[];
        Optden.AppendHistory=[];
        Optden.Slices=[];
        Optden.Frames=[];
        Optden.Overwrite=[];
        Optden.AdjustHeader='y';
        [~,~,~]=WriteBrik(brikout,HEAD,Optden);
        if Opts.brikout==2 
            ucom=['3dcopy ',filename,'_gsr_',alffout{loop1},'+tlrc.BRIK ' filename,'_gsr_',alffout{loop1},'.nii'];
            unix(ucom);
            ucomrm=['rm ',filename,'_gsr_',alffout{loop1},'+tlrc.*'];
            unix(ucomrm)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

