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
%   02/2010 AH: Replace unwrap_ifg_index with drop_ifg_index
%   01/2012 AH: Add bperp for new method 3D_NEW
%   01/2012 AH: Add back SULA error before unwrapping
%   02/2014 AH: Add predefined ph_uw option
%   10/2014 EH/DB: Suppress the aps removal correlated with topography
%   09/2015 DB: Include APS removal options
%   01/2017 DB: Double check if K_ps is not empty 
%   01/2017 DB: Check if the good values have right number of PS, if not
%               likely from an earlier pixel selection run, so do not apply.
%   06/2017 DB: Include stamps save for large variables
%   ======================================================================
logit;
fprintf('Phase-unwrapping...\n')
small_baseline_flag=getparm('small_baseline_flag',1);
unwrap_patch_phase=getparm('unwrap_patch_phase',1);
scla_deramp=getparm('scla_deramp',1);
subtr_tropo = getparm('subtr_tropo',1);
aps_name = getparm('tropo_method',1);


load psver
psname=['ps',num2str(psver)];
rcname=['rc',num2str(psver)];
pmname=['pm',num2str(psver)];
bpname=['bp',num2str(psver)];
goodname=['phuw_good',num2str(psver)];

if ~strcmpi(small_baseline_flag,'y')
    sclaname=['scla_smooth',num2str(psver)];
    apsname=['tca',num2str(psver)];
    phuwname=['phuw',num2str(psver),'.mat'];
else
    sclaname=['scla_smooth_sb',num2str(psver)];
    apsname=['tca_sb',num2str(psver)];
    phuwname=['phuw_sb',num2str(psver),'.mat'];
end

ps=load(psname);

drop_ifg_index=getparm('drop_ifg_index',1);
unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);

if exist(['./',bpname,'.mat'],'file')
    bp=load(bpname);
else
    bperp=ps.bperp;
    if ~strcmpi(small_baseline_flag,'y')
       bperp=bperp([1:ps.master_ix-1,ps.master_ix+1:end]);
    end
    bp.bperp_mat=repmat(bperp',ps.n_ps,1);
end

if ~strcmpi(small_baseline_flag,'y') 
    bperp_mat=[bp.bperp_mat(:,1:ps.master_ix-1),zeros(ps.n_ps,1,'single'),bp.bperp_mat(:,ps.master_ix:end)];
else
    bperp_mat=bp.bperp_mat;
end

if strcmpi(unwrap_patch_phase,'y')
    pm=load(pmname);
    ph_w=pm.ph_patch./abs(pm.ph_patch);
    clear pm
    if ~strcmpi(small_baseline_flag,'y')
        ph_w=[ph_w(:,1:ps.master_ix-1),ones(ps.n_ps,1),ph_w(:,ps.master_ix:end)];
    end
else
    rc=load(rcname);
    ph_w=rc.ph_rc;
    clear rc;
    if exist(['./',pmname,'.mat'],'file')
        pm=load(pmname,'K_ps');
        if isfield(pm,'K_ps')
            if ~isempty(pm.K_ps)
                ph_w=ph_w.*exp(j*(repmat(pm.K_ps,1,ps.n_ifg).*bperp_mat));
            end
        end
    end
end

ix=ph_w~=0;
ph_w(ix)=ph_w(ix)./abs(ph_w(ix)); % normalize, to avoid high freq artifacts being introduced in adaptive filtering

scla_subtracted_sw=0;
ramp_subtracted_sw=0;

options=struct('master_day',ps.master_day);
unwrap_hold_good_values=getparm('unwrap_hold_good_values',1);
if ~strcmpi(small_baseline_flag,'y') | ~exist(phuwname)
    unwrap_hold_good_values='n';
    logit('Code to hold good values skipped')
end
if unwrap_hold_good_values=='y'     
    sb_identify_good_pixels
    options.ph_uw_predef=nan(size(ph_w),'single');
    uw=load(phuwname);
    good=load(goodname);
    if ps.n_ps==size(good.good_pixels,1) & ps.n_ps==size(uw.ph_uw,1)
        options.ph_uw_predef(good.good_pixels)=uw.ph_uw(good.good_pixels);
    else
        fprintf('   wrong number of PS in keep good pixels - skipped...\n')
    end
    clear uw good;
end

if ~strcmpi(small_baseline_flag,'y') & exist([sclaname,'.mat'],'file') % PS
    fprintf('   subtracting scla and master aoe...\n')
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

if strcmpi(small_baseline_flag,'y') & exist([sclaname,'.mat'],'file') %Small baselines
    fprintf('   subtracting scla...\n')
    scla=load(sclaname);
    if size(scla.K_ps_uw,1)==ps.n_ps
      scla_subtracted_sw=1;
      ph_w=ph_w.*exp(-j*repmat(scla.K_ps_uw,1,ps.n_ifg).*bperp_mat); % subtract spatially correlated look angle error
      if unwrap_hold_good_values=='y'
          options.ph_uw_predef=options.ph_uw_predef-repmat(scla.K_ps_uw,1,ps.n_ifg).*bperp_mat; % subtract spatially correlated look angle error
      end
      if strcmpi(scla_deramp,'y') & isfield(scla,'ph_ramp') & size(scla.ph_ramp,1)==ps.n_ps
         ramp_subtracted_sw=1;
         ph_w=ph_w.*exp(-j*scla.ph_ramp); % subtract orbital ramps
         if unwrap_hold_good_values=='y'
             options.ph_uw_predef=options.ph_uw_predef-scla.ph_ramp;
         end
      end
    else
      fprintf('   wrong number of PS in scla - subtraction skipped...\n')
      delete([sclaname,'.mat'])
    end
    clear scla
end

clear bp

if exist([apsname,'.mat'],'file') && strcmpi(subtr_tropo,'y')
    fprintf('   subtracting slave aps...\n')
    aps=load(apsname);
    [aps_corr,fig_name_tca,aps_flag] = ps_plot_tca(aps,aps_name);

    ph_w=ph_w.*exp(-j*aps_corr);
    if unwrap_hold_good_values=='y'
        options.ph_uw_predef=options.ph_uw_predef-aps_corr;
    end
    clear aps
end


options.time_win=getparm('unwrap_time_win',1);
options.unwrap_method=getparm('unwrap_method',1);
options.grid_size=getparm('unwrap_grid_size',1);
options.prefilt_win=getparm('unwrap_gold_n_win',1);
options.goldfilt_flag=getparm('unwrap_prefilter_flag',1);
options.gold_alpha=getparm('unwrap_gold_alpha',1);
options.la_flag=getparm('unwrap_la_error_flag',1);
options.scf_flag=getparm('unwrap_spatial_cost_func_flag',1);



max_topo_err=getparm('max_topo_err',1);
lambda=getparm('lambda',1);


%%% ===============================================
%%% The code below needs to be made sensor specific
%%% ===============================================
rho = 830000; % mean range - need only be approximately correct
if isfield(ps,'mean_incidence')
    inc_mean=ps.mean_incidence;
else
    laname=['./la',num2str(psver),'.mat'];
    if exist(laname,'file')
        la=load(laname);
        inc_mean=mean(la.la)+0.052; % incidence angle approx equals look angle + 3 deg
        clear la
    else
        inc_mean=21*pi/180; % guess the incidence angle
    end
end
max_K=max_topo_err/(lambda*rho*sin(inc_mean)/4/pi);
%%% ===============================================
%%% The code above needs to be made sensor specific
%%% ===============================================



bperp_range=max(ps.bperp)-min(ps.bperp);
options.n_trial_wraps=(bperp_range*max_K/(2*pi));
logit(sprintf('n_trial_wraps=%f',options.n_trial_wraps))


if strcmpi(small_baseline_flag,'y')
    %options.lowfilt_flag='y';
    options.lowfilt_flag='n';
    ifgday_ix=ps.ifgday_ix;
    day=ps.day-ps.master_day;
else
    lowfilt_flag='n';
    %ifgday_ix=[];
    ifgday_ix=[ones(ps.n_ifg,1)*ps.master_ix,[1:ps.n_ifg]'];
    master_ix=sum(ps.master_day>ps.day)+1;
    unwrap_ifg_index=setdiff(unwrap_ifg_index,master_ix); % leave master ifg (which is only noise) out
    %day=ps.day(unwrap_ifg_index)-ps.master_day;
    day=ps.day-ps.master_day;
end

if unwrap_hold_good_values=='y'
    options.ph_uw_predef=options.ph_uw_predef(:,unwrap_ifg_index);
end

arch=computer('arch');
if ~strcmpi(arch(1:3),'win')
    [ph_uw_some,msd_some]=uw_3d(ph_w(:,unwrap_ifg_index),ps.xy,day,ifgday_ix(unwrap_ifg_index,:),ps.bperp(unwrap_ifg_index),options);
else
    logit('Windows detected: using old unwrapping code without statistical cost processing')
    [ph_uw_some]=uw_nosnaphu(ph_w(:,unwrap_ifg_index),ps.xy,day,options);
end


    
ph_uw=zeros(ps.n_ps,ps.n_ifg,'single');
msd=zeros(ps.n_ifg,1,'single');
ph_uw(:,unwrap_ifg_index)=ph_uw_some;
if exist('msd_some','var')
    msd(unwrap_ifg_index)=msd_some;
end

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

if exist([apsname,'.mat'],'file') && strcmpi(subtr_tropo,'y')
    fprintf('Adding back slave APS...\n')
    aps=load(apsname);
    [aps_corr,fig_name_tca,aps_flag] = ps_plot_tca(aps,aps_name);
    ph_uw=ph_uw+aps_corr;
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

stamps_save(phuwname,ph_uw,msd)
logit(1);
