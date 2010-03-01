function []=ps_sb_merge()
%PS_SB_MERGE merge PS and SB pixels
%
%   Andy Hooper May 2007
%


ps_wd=pwd;
i=strfind(ps_wd,'INSAR_');
if isempty(i)
    error('INSAR_* does not appear in current directory path')
end

sb_wd1=[ps_wd(1:i+13),'/SMALL_BASELINES'];
sb_wd=[sb_wd1,ps_wd(i+14:end)];
merged_wd1=[ps_wd(1:i+13),'/MERGED'];
merged_wd=[merged_wd1,ps_wd(i+14:end)];

if ~exist(merged_wd)
    mkdir(merged_wd)
end
if ~exist([merged_wd,'/patch_noover.in']) & exist([sb_wd,'/patch_noover.in'])
    aa=['!cp ',sb_wd,'/patch_noover.in ',merged_wd]; 
    eval(aa)
end
if ~exist([merged_wd1,'/parms.mat'])
    aa=['!cp ',sb_wd1,'/parms.mat ',merged_wd1,'/parms.mat']; 
    eval(aa)
end

psver=2;
save([merged_wd,'/psver'],'psver')

psps=load([ps_wd,'/ps2.mat']);
pssb=load([sb_wd,'/ps2.mat']);

pssb_ix=zeros(psps.n_ps,1); % points from ps to sb

imax=max([psps.ij(:,2);pssb.ij(:,2)])+1;
jmax=max([psps.ij(:,3);pssb.ij(:,3)])+1;
grid_ix=zeros(imax,jmax);
grid_ix(pssb.ij(:,2)+1+pssb.ij(:,3)*imax)=[1:pssb.n_ps];
pssb_ix=grid_ix(psps.ij(:,2)+1+psps.ij(:,3)*imax);


psu_ix=find(pssb_ix==0); % ps only index
psnu_ix=find(pssb_ix);  % ps non-unique index
sbnu_ix=pssb_ix(psnu_ix); % sb non-unique index

%%%%%%%%% Some non-adjacent pixels are allocated the same lon/lat by DORIS. %%%%%%%%%
[dups,I]=setdiff(psps.lonlat(psu_ix,:),pssb.lonlat,'rows');
fprintf('%d pixels with duplicate lon/lat dropped\n',length(psu_ix)-length(dups))
psu_ix=psu_ix(I);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    


G=zeros(pssb.n_ifg,pssb.n_image);
for i=1:pssb.n_ifg
     G(i,pssb.ifgday_ix(i,1))=-1;
     G(i,pssb.ifgday_ix(i,2))=1;
end

pmps=load([ps_wd,'/pm2.mat']);
pmsb=load([sb_wd,'/pm2.mat']);

pmps.ph_res_all=[pmps.ph_res(:,1:psps.master_ix-1),pmps.C_ps,pmps.ph_res(:,psps.master_ix:end)]; % add master residual
pmps.ph_res3=exp(j*(G*pmps.ph_res_all')'); % would-be residuals of sb ifgs
pmps.coh_ps2=abs(sum(pmps.ph_res3,2))/pssb.n_ifg; % would-be coherence of sb ifgs

psnu_coh=pmps.coh_ps2(psnu_ix);
psnu_snr=1./(1./psnu_coh.^2-1);
sbnu_coh=pmsb.coh_ps(sbnu_ix);
sbnu_snr=1./(1./sbnu_coh.^2-1);
ps_high_coh_ix=pmps.coh_ps2(psu_ix)>min(pmsb.coh_ps);
psu_ix=psu_ix(ps_high_coh_ix);
psu_coh=pmps.coh_ps2(psu_ix);


save([merged_wd,'/merge2'],'psu_ix','psu_coh','psnu_ix','sbnu_ix','psnu_coh','sbnu_coh','G')

ps=pssb;
ps.ij=[psps.ij(psu_ix,:);pssb.ij];
ps.lonlat=[psps.lonlat(psu_ix,:);pssb.lonlat];
ll0=(max(ps.lonlat)+min(ps.lonlat))/2;
xy=llh2local(ps.lonlat',ll0)*1000;
xy=xy';
sort_x=sortrows(xy,1);
sort_y=sortrows(xy,2);
n_pc=round(size(xy,1)*0.001);
bl=mean(sort_x(1:n_pc,:)); % bottom left corner
tr=mean(sort_x(end-n_pc:end,:)); % top right corner
br=mean(sort_y(1:n_pc,:)); % bottom right  corner
tl=mean(sort_y(end-n_pc:end,:)); % top left corner

heading=getparm('heading');
theta=(180-heading)*pi/180;
if theta>pi
    theta=theta-2*pi;
end

rotm=[cos(theta),sin(theta); -sin(theta),cos(theta)];
xy=xy';
xynew=rotm*xy; % rotate so that scene axes approx align with x=0 and y=0
if max(xynew(1,:))-min(xynew(1,:))<max(xy(1,:))-min(xy(1,:)) &...
   max(xynew(2,:))-min(xynew(2,:))<max(xy(2,:))-min(xy(2,:))
    xy=xynew; % check that rotation is an improvement
    disp(['Rotating by ',num2str(theta*180/pi),' degrees']);
end

xy=single(xy');
[dummy,sort_ix]=sortrows(xy,[2,1]); % sort in ascending y order
xy=xy(sort_ix,:);
xy=[[1:size(xy,1)]',xy];
xy(:,2:3)=round(xy(:,2:3)*1000)/1000; % round to mm
ps.xy=xy;
ps.lonlat=ps.lonlat(sort_ix,:);
ps.ij=ps.ij(sort_ix,:);
ps.n_ps=size(ps.xy,1);
save([merged_wd,'/ps2'],'-struct','ps')

pmps.ph_patch2=exp(j*G*angle([pmps.ph_patch(:,1:psps.master_ix-1),ones(psps.n_ps,1),pmps.ph_patch(:,psps.master_ix:end)])').';
pm.ph_patch=pmsb.ph_patch;
pm.ph_patch(sbnu_ix,:)=pmsb.ph_patch(sbnu_ix,:).*repmat(sbnu_snr,1,pssb.n_ifg)+pmps.ph_patch2(psnu_ix,:).*repmat(psnu_snr,1,pssb.n_ifg);
pm.ph_patch(sbnu_ix,:)=pm.ph_patch(sbnu_ix,:)./abs(pm.ph_patch(sbnu_ix,:));
pm.ph_patch=[pmps.ph_patch2(psu_ix,:);pm.ph_patch];
pm.ph_patch=pm.ph_patch(sort_ix,:);

save([merged_wd,'/pm2'],'-struct','pm')
clear pm pmps pmsb

rcps=load([ps_wd,'/rc2.mat']);
rcsb=load([sb_wd,'/rc2.mat']);
rcps.ph_rc2=exp(j*G*angle(rcps.ph_rc)').'; % convert phase to sb ifgs
rc=rcsb;
rc.ph_rc(sbnu_ix,:)=rcsb.ph_rc(sbnu_ix,:).*repmat(sbnu_snr,1,pssb.n_ifg)+rcps.ph_rc2(psnu_ix,:).*repmat(psnu_snr,1,pssb.n_ifg);
rc.ph_rc(sbnu_ix,:)=rc.ph_rc(sbnu_ix,:)./abs(rc.ph_rc(sbnu_ix,:));
%rc.ph_rc(sbnu_ix,:)=rc.ph_rc(sbnu_ix,:)+rcps.ph_rc2(psnu_ix,:);
rc.ph_rc=[rcps.ph_rc2(psu_ix,:);rc.ph_rc];
rc.ph_rc=rc.ph_rc(sort_ix,:);

save([merged_wd,'/rc2'],'-struct','rc')
clear rc rcps rcsb


bpps=load([ps_wd,'/bp2.mat']);
bpsb=load([sb_wd,'/bp2.mat']);
bpps.bperp_mat2=(G*[bpps.bperp_mat(:,1:psps.master_ix-1),zeros(psps.n_ps,1),bpps.bperp_mat(:,psps.master_ix:end)]')';
bp=bpsb;
bp.bperp_mat=[bpps.bperp_mat2(psu_ix,:);bp.bperp_mat];
bp.bperp_mat=bp.bperp_mat(sort_ix,:);
save([merged_wd,'/bp2'],'-struct','bp')
clear bp bpps bpsb


hgtps=load([ps_wd,'/hgt2.mat']);
hgtsb=load([sb_wd,'/hgt2.mat']);
hgt=hgtsb;
hgt.hgt=[hgtps.hgt(psu_ix,:);hgt.hgt];
hgt.hgt=hgt.hgt(sort_ix,:);
save([merged_wd,'/hgt2'],'-struct','hgt')
clear hgt hgtps hgtsb

laps=load([ps_wd,'/la2.mat']);
lasb=load([sb_wd,'/la2.mat']);
la=lasb;
la.la=[laps.la(psu_ix,:);la.la];
la.la=la.la(sort_ix,:);
save([merged_wd,'/la2'],'-struct','la')
clear la laps lasb

