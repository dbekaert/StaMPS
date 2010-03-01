function []=sb_invert_uw()
%SB_INVERT_UW Invert unwrapped phase of short baseline ifgs
%
%   Andy Hooper, September 2006



load psver
psname=['ps',num2str(psver)];
phuwsbname=['phuw_sb',num2str(psver)];
phuwname=['phuw',num2str(psver)];

phuwsb=load(phuwsbname);

unwrap_ifg_index=getparm('unwrap_ifg_index');
if strcmpi(unwrap_ifg_index,'all')
    unwrap_ifg_index=[1:ps.n_ifg];
end

ph_uw_sb=phuwsb.ph_uw(:,unwrap_ifg_index);
clear phuwsb

ps=load(psname);
ifgday_ix=ps.ifgday_ix(unwrap_ifg_index,:);
n_ifg=length(unwrap_ifg_index);

G=zeros(n_ifg,ps.n_image);
for i=1:n_ifg
    G(i,ifgday_ix(i,1))=-1;
    G(i,ifgday_ix(i,2))=1;
end
G=G(:,[1:ps.master_ix-1,ps.master_ix+1:end]); % take out master as ref

ph_uw=zeros(ps.n_ps,ps.n_image-1);
for i=1:ps.n_ps
    ph_uw(i,:)=(G\ph_uw_sb(i,:)')';
end

ph_uw=[ph_uw(:,1:ps.master_ix-1),zeros(ps.n_ps,1),ph_uw(:,ps.master_ix:end)];
save(phuwname,'ph_uw')

