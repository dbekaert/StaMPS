function []=sb_invert_uw()
%SB_INVERT_UW Invert unwrapped phase of short baseline ifgs
%
%   Andy Hooper, September 2006
%
%   ======================================================================
%   07/2008 AH: Allow isolated images, not included in any interferogram
%   11/2009 AH: Weight inversion using var-covar matrix
%   01/2010 AH: Correct lscov weighting and ensure well-conditioned
%   02/2010 AH: Reference phase to reference area before inverting
%   02/2010 AH: Replace unwrap_ifg_index with drop_ifg_index
%   03/2010 AH: Save small baseline covariance
%   07/2010 AH: Default small baseline covariance for multilooked case
%   01/2013 AH: Invert topo-correlated atmosphere if calculated
%   04/2013 AH: Remove inversion of topo-correlated atmosphere
%   09/2014 DB: Fix for nanmean of the reference area
%   06/2017 DB: Include stamps save for large variables
%   ======================================================================
logit;

load psver
psname=['./ps',num2str(psver)];
rcname=['./rc',num2str(psver)];
pmname=['./pm',num2str(psver)];
phuwsbname=['./phuw_sb',num2str(psver)];
phuwsbresname=['./phuw_sb_res',num2str(psver)];
phuwname=['./phuw',num2str(psver)];
stratsbname=['./tca_sb',num2str(psver)];
stratname=['./tca',num2str(psver)];

ps=load(psname);

drop_ifg_index=getparm('drop_ifg_index');
unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);

if exist([pmname,'.mat'],'file') 
    rc=load(rcname);
    pm=load(pmname);
    
    if isfield(pm,'ph_patch')
        if ~isempty(pm.ph_patch)
            ph_noise=angle(rc.ph_rc.*conj(pm.ph_patch));
            clear pm rc
            sb_cov=double(cov(ph_noise)); % Covariance between IFGs
            clear ph_noise
        else
            sb_cov=eye(ps.n_ifg);
        end
    else
        sb_cov=eye(ps.n_ifg);
    end
else
    sb_cov=eye(ps.n_ifg);
end
C=sb_cov(unwrap_ifg_index,unwrap_ifg_index);

phuwsb=load(phuwsbname);
ph_uw_sb=phuwsb.ph_uw(:,unwrap_ifg_index);
clear phuwsb

ref_ps=ps_setref;

%if  exist([stratsbname,'.mat'],'file')
%    strat=load(stratsbname);
%    if size(strat.strat_corr,1)==ps.n_ps
%        strat_sb=strat.strat_corr(:,unwrap_ifg_index);
%        clear strat
%        strat_sb=strat_sb-repmat(mean(strat_sb(ref_ps,:)),ps.n_ps,1);
%    end
%end

ph_uw_sb=ph_uw_sb-repmat(nanmean(ph_uw_sb(ref_ps,:)),ps.n_ps,1);

G=zeros(ps.n_ifg,ps.n_image);
for i=1:ps.n_ifg
    G(i,ps.ifgday_ix(i,1))=-1;
    G(i,ps.ifgday_ix(i,2))=1;
end
%G=G(:,[1:ps.master_ix-1,ps.master_ix+1:end]); % take out master as ref
if sum(abs(G(:,ps.master_ix)))==0; 
   error('Apparently none of the unwrapped interferograms include the original master image')
else
    G(:,ps.master_ix)=0; % take out master as ref by setting to zero
end


G2=G(unwrap_ifg_index,:);

nzc_ix=sum(abs(G2))~=0; % index for non-zero columns
G2=G2(:,nzc_ix);


if rank(G2)<size(G2,2) 
    stamps_save(phuwsbresname,sb_cov)
    error('There are isolated subsets (cannot be inverted w.r.t. master)')
end

ph_uw=zeros(ps.n_ps,ps.n_image,'single');
%ph_uw(:,nzc_ix)=(G2\double(ph_uw_sb'))';
while rcond(C)<0.001 % ensure not close to singular
    C=C+eye(size(C,1))*0.01;
end
ph_uw(:,nzc_ix)=lscov(G2,double(ph_uw_sb'),C)';
clear ph_uw_sm
sm_cov=zeros(ps.n_image);
sm_cov(nzc_ix,nzc_ix)=inv(G2'*inv(C)*G2);

ph_res=single(G*ph_uw')';

%ph_uw=[ph_uw(:,1:ps.master_ix-1),zeros(ps.n_ps,1),ph_uw(:,ps.master_ix:end)];

unwrap_ifg_index_sm=[1:ps.n_image]; % single master index
nzc_ix(ps.master_ix)=1;
unwrap_ifg_index_sm=unwrap_ifg_index_sm(nzc_ix);

stamps_save(phuwname,ph_uw,unwrap_ifg_index_sm)
stamps_save(phuwsbresname,ph_res,sb_cov,sm_cov)

%if exist('strat_sb','var')
%    strat_corr=zeros(size(ph_uw),'single');
%    strat_corr(:,nzc_ix)=lscov(G2,double(strat_sb'))';
%    save(stratname,'strat')
%end

logit(1);
