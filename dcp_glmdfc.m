function [er_betas,er_tvals,er_pvals]=dcp_glmdfc(indvs,dmat)

% Code by Michael J. Tobia, Ph.D. as part of the 
% Dynamic Connectivity Processing (DCP) toolbox
% DCP_v1.01 private release 1/10/2018

% indvs is time dependent 3D tensor from a single subject
% dmat is predictor variables with same length as time, i.e., a design matrix
% outputs are 3D tensors of betas, tvals, and pvals

outtens=0;

[~,rr,~]=size(indvs);
OM=dcp_ten2mat(indvs);
[~,c1]=size(OM);
[~,c2]=size(dmat);
betas1=zeros(c1,c2+1);
pvals1=zeros(c1,c2+1);
tvals1=zeros(c1,c2+1);

wb=waitbar(0,'Calculating Multiple Regression...');
for loop1=1:c1
    waitbar(loop1/c1,wb,'Calculating Multiple Regression...');
    statsout=regstats(OM(:,loop1),dmat,'linear','tstat');
    betas1(loop1,:)=statsout.tstat.beta';
    pvals1(loop1,:)=statsout.tstat.pval';
    tvals1(loop1,:)=statsout.tstat.t';
end

if outtens==1
    er_betas=zeros(c2+1,rr,rr);
    er_tvals=zeros(c2+1,rr,rr);
    er_pvals=zeros(c2+1,rr,rr);
    for loop2=1:c2+1
        waitbar(loop2/c2+1,wb,'Calculating Multiple Regression...');
        er_betas(loop2,:,:)=dcp_mat2tens(betas1(:,loop2));
        er_tvals(loop2,:,:)=dcp_mat2tens(tvals1(:,loop2));
        er_pvals(loop2,:,:)=dcp_mat2tens(pvals1(:,loop2));    
    end
else
    er_betas=betas1;
    er_tvals=tvals1;
    er_pvals=pvals1;    
end
waitbar(1/1,wb,'Done!!');
close(wb)

end

