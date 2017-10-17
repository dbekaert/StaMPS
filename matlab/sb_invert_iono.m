function []=sb_invert_iono(iono_flag)
%SB_invert_iono Invert unwrapped phase of short baseline ifgs
%
%   Bekaert David, November 2014   - University of Leeds
%
% modifications:
% 12/2014   DB  Invert such that when one date exist the delay is excluded

logit;

load psver
psname=['./ps',num2str(psver)];
iononame = ['./ica' num2str(psver) '.mat'];
ionosbname = ['./ica_sb' num2str(psver) '.mat'];


if nargin<1
   iono_flag = []; 
end
ps=load(psname);

drop_ifg_index=getparm('drop_ifg_index');
unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);

ref_ps=ps_setref;

G=zeros(ps.n_ifg,ps.n_image);
for i=1:ps.n_ifg
    G(i,ps.ifgday_ix(i,1))=-1;
    G(i,ps.ifgday_ix(i,2))=1;
end
if sum(abs(G(:,ps.master_ix)))==0; 
   error('Apparently none of the unwrapped interferograms include the original master image')
else
    G(:,ps.master_ix)=0; % take out master as ref by setting to zero
end


G2=G(unwrap_ifg_index,:);
nzc_ix=sum(abs(G2))~=0; % index for non-zero columns
G2=G2(:,nzc_ix);
if rank(G2)<size(G2,2) 
    error('There are isolated subsets (cannot be inverted w.r.t. master)')
end


if ~isempty(iono_flag)
     aps=load(ionosbname);
    [aps_corr_sb,fig_name_tca] = ps_plot_ica(aps,iono_flag);
    
    
    ix = find(nanmean(aps_corr_sb,1)==0);
    unwrap_ifg_index_new = setdiff(unwrap_ifg_index,ix);

    G3=G(unwrap_ifg_index_new,:);
    nzc_ix=sum(abs(G3))~=0; % index for non-zero columns
    G3=G3(:,nzc_ix);
    if rank(G3)<size(G3,2) 
        error('There are isolated subsets that do not have an APS estimate (cannot be inverted w.r.t. master)')
    end

    aps_corr_sb= aps_corr_sb(:,unwrap_ifg_index_new);
    G2 =G3;
    
    % check if the master to be inverted actually has a delay estimated with it.
    ifgday_new = ps.ifgday(unwrap_ifg_index_new,:);
    ifgs_days = unique(reshape(ifgday_new((sum(aps_corr_sb)~=0),:),[],1));

    if sum(ps.master_day==ifgs_days)~=1
        error('Master does not have an APS estimated, inversion not possible')
    end

    aps_corr=zeros(ps.n_ps,ps.n_image,'single');
    aps_corr(:,nzc_ix)=lscov(G2,double(aps_corr_sb'))';
    
    
    
    if iono_flag==1 % linear correction
        ph_iono_azshift = aps_corr;       
        if exist(iononame,'file')==2
            save(iononame,'-append','ph_iono_azshift')
        else
            save(iononame,'ph_iono_azshift')
        end  

    end
end



logit(1);
