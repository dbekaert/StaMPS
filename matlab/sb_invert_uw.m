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
%   ======================================================================


load psver
psname=['ps',num2str(psver)];
rcname=['rc',num2str(psver)];
pmname=['pm',num2str(psver)];
phuwsbname=['phuw_sb',num2str(psver)];
phuwsbresname=['phuw_sb_res',num2str(psver)];
phuwname=['phuw',num2str(psver)];

ps=load(psname);

unwrap_ifg_index=getparm('unwrap_ifg_index');
if strcmpi(unwrap_ifg_index,'all')
    unwrap_ifg_index=[1:ps.n_ifg];
end

rc=load(rcname);
pm=load(pmname);
ph_noise=angle(rc.ph_rc.*conj(pm.ph_patch));
clear pm rc
ph_noise=ph_noise(:,unwrap_ifg_index);
C=double(cov(ph_noise)); % Covariance between IFGs
clear ph_noise

phuwsb=load(phuwsbname);
ph_uw_sb=phuwsb.ph_uw(:,unwrap_ifg_index);
clear phuwsb

ref_ps=ps_setref;
ph_uw_sb=ph_uw_sb-repmat(mean(ph_uw_sb(ref_ps,:)),ps.n_ps,1);

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
   error('There are isolated images (cannot be inverted w.r.t. master)')
end

ph_uw=zeros(ps.n_ps,ps.n_image,'single');
%ph_uw(:,nzc_ix)=(G2\double(ph_uw_sb'))';
while rcond(C)<0.001 % ensure not close to singular
    C=C+eye(size(C,1))*0.01;
end
ph_uw(:,nzc_ix)=lscov(G2,double(ph_uw_sb'),C)';

ph_res=single(G*ph_uw')';

%ph_uw=[ph_uw(:,1:ps.master_ix-1),zeros(ps.n_ps,1),ph_uw(:,ps.master_ix:end)];

unwrap_ifg_index_sm=[1:ps.n_image]; % single master index
nzc_ix(ps.master_ix)=1;
unwrap_ifg_index_sm=unwrap_ifg_index_sm(nzc_ix);

save(phuwname,'ph_uw','unwrap_ifg_index_sm')
save(phuwsbresname,'ph_res')

