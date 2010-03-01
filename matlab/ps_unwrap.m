function []=ps_unwrap()
%PS_UNWRAP unwrap phase using the 3-D cost function phase unwrapping algorithm 
%
%   Andy Hooper, Jun 2006
%
%   ======================================================================
%   03/2008 AH: Fixed bug to make prefilter='n' option work
%   03/2009 AH: Use smoothed scla instead of scla files
%   06/2009 AH: Orbital ramps option added 
%   08/2009 AH: Goldstein alpha value added to options
%   ======================================================================

fprintf('Phase-unwrapping...\n')

small_baseline_flag=getparm('small_baseline_flag',1);
unwrap_patch_phase=getparm('unwrap_patch_phase',1);
scla_deramp=getparm('scla_deramp',1);

load psver
psname=['ps',num2str(psver)];
rcname=['rc',num2str(psver)];
pmname=['pm',num2str(psver)];
bpname=['bp',num2str(psver)];
if ~strcmpi(small_baseline_flag,'y')
    sclaname=['scla_smooth',num2str(psver)];
    apsname=['aps',num2str(psver)];
    phuwname=['phuw',num2str(psver)];
else
    sclaname=['scla_smooth_sb',num2str(psver)];
    apsname=['aps_sb',num2str(psver)];
    phuwname=['phuw_sb',num2str(psver)];
end

ps=load(psname);

unwrap_ifg_index=getparm('unwrap_ifg_index');
unwrap_ifg_index=sort(unwrap_ifg_index);
if strcmp(unwrap_ifg_index,'all')
    unwrap_ifg_index=[1:ps.n_ifg];
end

if strcmpi(unwrap_patch_phase,'y')
    pm=load(pmname);
    ph_w=pm.ph_patch./abs(pm.ph_patch);
    clear pm
    if ~strcmpi(small_baseline_flag,'y')
        ph_w=[ph_w(:,1:ps.master_ix-1),zeros(ps.n_ps,1),ph_w(:,ps.master_ix:end)];
    end
else
    rc=load(rcname);
    ph_w=rc.ph_rc;
    clear rc;
end

ix=ph_w~=0;
ph_w(ix)=ph_w(ix)./abs(ph_w(ix)); % normalize, to avoid high freq artifacts being introduced in adaptive filtering


if exist(['./',bpname,'.mat'],'file')
    bp=load(bpname);
else
    bperp=ps.bperp;
    if ~strcmpi(small_baseline_flag,'y')
       bperp=bperp([1:ps.master_ix-1,ps.master_ix+1:end]);
    end
    bp.bperp_mat=repmat(bperp',ps.n_ps,1);
end

scla_subtracted_sw=0;
ramp_subtracted_sw=0;

if ~strcmpi(small_baseline_flag,'y') & exist([sclaname,'.mat'],'file')
    fprintf('   subtracting scla and master aoe...\n')
    bperp_mat=[bp.bperp_mat(:,1:ps.master_ix-1),zeros(ps.n_ps,1,'single'),bp.bperp_mat(:,ps.master_ix:end)];
    scla=load(sclaname);
    if size(scla.K_ps_uw,1)==ps.n_ps
      scla_subtracted_sw=1;
      ph_w=ph_w.*exp(-j*repmat(scla.K_ps_uw,1,ps.n_ifg).*bperp_mat); % subtract spatially correlated look angle error
      ph_w=ph_w.*repmat(exp(-j*scla.C_ps_uw),1,ps.n_ifg); % subtract master APS
      if strcmpi(scla_deramp,'y') & isfield(scla,'ph_ramp') & size(scla.ph_ramp,1)==ps.n_ps
         ramp_subtracted_sw=1;
         ph_w=ph_w.*exp(-j*scla.ph_ramp); % subtract orbital ramps
      end
    else
      fprintf('   wrong number of PS in scla - subtraction skipped...\n')
      delete([sclaname,'.mat'])
    end
    clear scla
end

if strcmpi(small_baseline_flag,'y') & exist([sclaname,'.mat'],'file')
    fprintf('   subtracting scla...\n')
    bperp_mat=bp.bperp_mat;
    scla=load(sclaname);
    if size(scla.K_ps_uw,1)==ps.n_ps
      scla_subtracted_sw=1;
      ph_w=ph_w.*exp(-j*repmat(scla.K_ps_uw,1,ps.n_ifg).*bperp_mat); % subtract spatially correlated look angle error
      if strcmpi(scla_deramp,'y') & isfield(scla,'ph_ramp') & size(scla.ph_ramp,1)==ps.n_ps
         ramp_subtracted_sw=1;
         ph_w=ph_w.*exp(-j*scla.ph_ramp); % subtract orbital ramps
      end
    else
      fprintf('   wrong number of PS in scla - subtraction skipped...\n')
      delete([sclaname,'.mat'])
    end
    clear scla
end

clear bp

if exist([apsname,'.mat'],'file')
    fprintf('   subtracting slave aps...\n')
    aps=load(apsname);
    ph_w=ph_w.*exp(-j*aps.ph_aps_slave);
    clear aps
end


options=struct('master_day',ps.master_day);
options.time_win=getparm('unwrap_time_win',1);
options.unwrap_method=getparm('unwrap_method',1);
options.grid_size=getparm('unwrap_grid_size',1);
options.prefilt_win=getparm('unwrap_gold_n_win',1);
options.goldfilt_flag=getparm('unwrap_prefilter_flag',1);
options.gold_alpha=getparm('unwrap_gold_alpha',1);

if strcmpi(small_baseline_flag,'y')
    %options.lowfilt_flag='y';
    options.lowfilt_flag='n';
    ifgday_ix=ps.ifgday_ix;
    day=ps.day-ps.master_day;
else
    lowfilt_flag='n';
    ifgday_ix=[];
    master_ix=sum(ps.master_day>ps.day)+1;
    unwrap_ifg_index=setdiff(unwrap_ifg_index,master_ix); % leave master ifg (which is only noise) out
    day=ps.day(unwrap_ifg_index)-ps.master_day;
end

[ph_uw_some,msd_some]=uw_3d(ph_w(:,unwrap_ifg_index),ps.xy,day,ifgday_ix,options);

ph_uw=zeros(ps.n_ps,ps.n_ifg,'single');
msd=zeros(ps.n_ifg,1,'single');
ph_uw(:,unwrap_ifg_index)=ph_uw_some;
msd(unwrap_ifg_index)=msd_some;

if scla_subtracted_sw & ~strcmpi(small_baseline_flag,'y')
    fprintf('Adding back SCLA and master AOE...\n')
    scla=load(sclaname);
    ph_uw=ph_uw+(repmat(scla.K_ps_uw,1,ps.n_ifg).*bperp_mat); % add back spatially correlated look angle error
    ph_uw=ph_uw+repmat(scla.C_ps_uw,1,ps.n_ifg); % add back master APS
    if ramp_subtracted_sw
        ph_uw=ph_uw+scla.ph_ramp; % add back orbital ramps
    end
    clear bp scla
end

if scla_subtracted_sw & strcmpi(small_baseline_flag,'y')
    fprintf('Adding back SCLA...\n')
    scla=load(sclaname);
    ph_uw=ph_uw+(repmat(scla.K_ps_uw,1,ps.n_ifg).*bperp_mat); % add back spatially correlated look angle error
    if ramp_subtracted_sw
        ph_uw=ph_uw+scla.ph_ramp; % add back orbital ramps
    end
    clear bp scla
end

if exist([apsname,'.mat'],'file')
    fprintf('Adding back slave APS...\n')
    aps=load(apsname);
    ph_uw=ph_uw+aps.ph_aps_slave;
    clear aps
end

if strcmpi(unwrap_patch_phase,'y')
    pm=load(pmname);
    ph_w=pm.ph_patch./abs(pm.ph_patch);
    clear pm
    if ~strcmpi(small_baseline_flag,'y')
        ph_w=[ph_w(:,1:ps.master_ix-1),zeros(ps.n_ps,1),ph_w(:,ps.master_ix:end)];
    end
    rc=load(rcname);
    ph_uw=ph_uw+angle(rc.ph_rc.*conj(ph_w));
end

ph_uw(:,setdiff([1:ps.n_ifg],unwrap_ifg_index))=0;

save(phuwname,'ph_uw','msd')
