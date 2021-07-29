function []=ps_calc_scla(use_small_baselines,coest_mean_vel)
% PS_CALC_SCLA calculate SCLA and master atmosphere & orbit error
%    PS_CALC_SCLA(USE_SMALL_BASELINES) Set USE_SMALL_BASELINES to 1 
%    to use small baseline interferograms
%
%   Andy Hooper, Nov 2006
%
%   ================================================================
%   01/2008 AH: Date processing changed for non-english locales
%   06/2009 AH: Orbital ramps option added 
%   06/2009 AH: Ramps option for small baselines corrected 
%   02/2010 AH: Replace unwrap_ifg_index with drop_ifg_index
%   03/2010 AH: Include var/cov in inversion
%   05/2010 AH: add ./ before .mat files 
%   05/2010 AH: check for ifgstd and don't solve for vel if < 4 ifgs 
%   05/2010 AH: correct ramp processing fro dropped ifgs
%   08/2010 AH: use recalc_index correctly in small baseline case
%   08/2010 AH: use mean bperp value
%   11/2010 AH: replace recalc_index with scla_drop_index
%   01/2012 AH: use short time sep ifgs for single master SCLA calc
%   07/2014 AH: remove Rank Deficient warning 
%   09/2015 DB: Include TRAIN support
%   09/2015 DB/EH: Debug nans, deramping fix using script ps_deramp
%   09/2016 AH: Drop master from single master baseline calc for SB
%   06/2017 DB: include stamps_save for larger variables
%   09/2017 DB: drop_master_ix is not defined, fix this
%   ================================================================
logit;
logit(sprintf('Estimating spatially-correlated look angle error...'),2)

if nargin<1
    use_small_baselines=0;
end
if nargin<2
    coest_mean_vel=0;
end


small_baseline_flag=getparm('small_baseline_flag',1);
drop_ifg_index=getparm('drop_ifg_index',1);
scla_method=getparm('scla_method',1);
scla_deramp=getparm('scla_deramp',1);
subtr_tropo=getparm('subtr_tropo',1);
tropo_method=getparm('tropo_method',1);


if use_small_baselines~=0
    if small_baseline_flag~='y'
        error('   Use small baselines requested but there are none')
    end
end

if use_small_baselines==0
    scla_drop_index=getparm('scla_drop_index',1);
else
    scla_drop_index=getparm('sb_scla_drop_index',1);
    fprintf('   Using small baseline interferograms\n')
end

load psver
psname=['./ps',num2str(psver)];
rcname=['./rc',num2str(psver)];
pmname=['./pm',num2str(psver)];
bpname=['./bp',num2str(psver)];
meanvname=['./mv',num2str(psver)];
ifgstdname=['./ifgstd',num2str(psver)];
phuwsbresname=['./phuw_sb_res',num2str(psver)];
if use_small_baselines==0
    phuwname=['./phuw',num2str(psver)];
    sclaname=['./scla',num2str(psver)];
    apsname_old=['./aps',num2str(psver)];           % renamed to old
    apsname=['./tca',num2str(psver)];               % the new tca option
else
    phuwname=['./phuw_sb',num2str(psver)];
    sclaname=['./scla_sb',num2str(psver)];
    apsname_old=['./aps_sb',num2str(psver)];        % renamed to old
    apsname=['./tca_sb',num2str(psver)];            % the new tca option
end


if use_small_baselines==0
    evalcmd=['!rm -f ',meanvname,'.mat'];
    eval(evalcmd)
end

ps=load(psname);
if exist([bpname,'.mat'],'file')
    bp=load(bpname);
else
    bperp=ps.bperp;
    if ~strcmpi(small_baseline_flag,'y')
       bperp=bperp([1:ps.master_ix-1,ps.master_ix+1:end]);
    end
    bp.bperp_mat=repmat(bperp',ps.n_ps,1);
end
uw=load(phuwname);

if strcmpi(small_baseline_flag,'y') & use_small_baselines==0
    unwrap_ifg_index=[1:ps.n_image];
    n_ifg=ps.n_image;
else
    unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);
    n_ifg=ps.n_ifg;
end


if strcmpi(subtr_tropo,'y')
    % Remove the tropo correction - TRAIN support
    % recompute the APS inversion on the fly as user migth have dropped
    % SB ifgs before and needs new update of the SM APS too.
    
%    if exist(apsname,'file')~=2
%        % the tca file does not exist. See in case this is SM if it needs
%        % to be inverted 
        if strcmpi(apsname,['./tca',num2str(psver)])
            if strcmpi(getparm('small_baseline_flag'),'y')
                sb_invert_aps(tropo_method);
            end
        end
%    end
    aps = load(apsname);
    [aps_corr,fig_name_tca,tropo_method] = ps_plot_tca(aps,tropo_method);
    uw.ph_uw=uw.ph_uw-aps_corr;
end

if strcmpi(scla_deramp,'y')
    fprintf('\n   deramping ifgs...\n')
    
    [ph_all,ph_ramp] = ps_deramp(ps,uw.ph_uw);
    uw.ph_uw = uw.ph_uw - ph_ramp;
    
%     ph_ramp=zeros(ps.n_ps,n_ifg,'single');
%     G=double([ones(ps.n_ps,1),ps.xy(:,2),ps.xy(:,3)]);
%     for i=1:length(unwrap_ifg_index)
%         d=uw.ph_uw(:,unwrap_ifg_index(i));
%         m=G\double(d(:));
%         ph_this_ramp=G*m;
%         uw.ph_uw(:,unwrap_ifg_index(i))=uw.ph_uw(:,unwrap_ifg_index(i))-ph_this_ramp; % subtract ramp
%         ph_ramp(:,unwrap_ifg_index(i))=ph_this_ramp;
%     end
else
    ph_ramp=[];
end

unwrap_ifg_index=setdiff(unwrap_ifg_index,scla_drop_index);


% Check with Andy:
% 1) should this not be placed before the ramp computation.
% 2) if this is spatial fitlering in time - not compatible with TRAIN
if exist([apsname_old,'.mat'],'file')
    if strcmpi(subtr_tropo,'y')
        fprintf(['You are removing atmosphere twice. Do not do this, either do:\n use ' apsname_old ' with subtr_tropo=''n''\n remove ' apsname_old ' use subtr_tropo=''y''\n'])
    end
    aps=load(apsname_old);
    uw.ph_uw=uw.ph_uw-aps.ph_aps_slave;
end

ref_ps=ps_setref;
uw.ph_uw=uw.ph_uw-repmat(nanmean(uw.ph_uw(ref_ps,:),1),ps.n_ps,1);

if use_small_baselines==0
    if strcmpi(small_baseline_flag,'y')
        bperp_mat=zeros(ps.n_ps,ps.n_image,'single');
        G=zeros(ps.n_ifg,ps.n_image);
        for i=1:ps.n_ifg
             G(i,ps.ifgday_ix(i,1))=-1;
             G(i,ps.ifgday_ix(i,2))=1;
        end
        if isfield(uw,'unwrap_ifg_index_sm')
            unwrap_ifg_index=setdiff(uw.unwrap_ifg_index_sm,scla_drop_index);
        end
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);

        G=G(:,unwrap_ifg_index);
        bperp_some=[G\double(bp.bperp_mat')]';
        bperp_mat(:,unwrap_ifg_index)=bperp_some;
        clear bperp_some
    else
        bperp_mat=[bp.bperp_mat(:,1:ps.master_ix-1),zeros(ps.n_ps,1,'single'),bp.bperp_mat(:,ps.master_ix:end)];
    end
    day=diff(ps.day(unwrap_ifg_index));
    ph=double(diff(uw.ph_uw(:,unwrap_ifg_index),[],2)); % sequential dph, to reduce influence of defo
    bperp=diff(bperp_mat(:,unwrap_ifg_index),[],2);
else
    
    bperp_mat=bp.bperp_mat;
    bperp=bperp_mat(:,unwrap_ifg_index);
    day=ps.ifgday(unwrap_ifg_index,2)-ps.ifgday(unwrap_ifg_index,1);
    ph=double(uw.ph_uw(:,unwrap_ifg_index));
end
clear bp

bprint=mean(bperp);
logit(sprintf('%d ifgs used in estimation:',size(ph,2)))

for i=1:size(ph,2)
    if use_small_baselines~=0
        logit(sprintf('   %s to %s %5d days %5d m',datestr(ps.ifgday(unwrap_ifg_index(i),1)),datestr(ps.ifgday(unwrap_ifg_index(i),2)),day(i),round(bprint(i))))
    else
        logit(sprintf('   %s to %s %5d days %5d m',datestr(ps.day(unwrap_ifg_index(i))),datestr(ps.day(unwrap_ifg_index(i+1))),day(i),round(bprint(i))))
    end
end

K_ps_uw=zeros(ps.n_ps,1);

if coest_mean_vel==0 | length(unwrap_ifg_index)<4
    G=[ones(size(ph,2),1),double(mean(bperp)')];
else
    G=[ones(size(ph,2),1),double(mean(bperp)'),double(day)];
end

ifg_vcm=eye(ps.n_ifg);
    
if strcmpi(small_baseline_flag,'y')
    if use_small_baselines==0 
        phuwres=load(phuwsbresname,'sm_cov');
        if isfield(phuwres,'sm_cov')
            ifg_vcm=phuwres.sm_cov;
        end
    else
        phuwres=load(phuwsbresname,'sb_cov');
        if isfield(phuwres,'sb_cov')
            ifg_vcm=phuwres.sb_cov;
        end
    end
else
    if exist([ifgstdname,'.mat'])
        ifgstd=load(ifgstdname);
        ifg_vcm=double(diag((ifgstd.ifg_std*pi/180).^2));
        clear ifgstd
    end
end

if use_small_baselines==0 
    ifg_vcm_use=eye(size(ph,2)); % don't know true var/cov because of non-lin motion and APS
else
    ifg_vcm_use=ifg_vcm(unwrap_ifg_index,unwrap_ifg_index);
end


m=lscov(G,ph',ifg_vcm_use); % L2-norm
K_ps_uw=m(2,:)';
if coest_mean_vel~=0
    v_ps_uw=m(3,:)';
end

if strcmpi(scla_method,'L1')
    for i=1:ps.n_ps
        d=ph(i,:)';
        m2=m(:,i);
        m2=fminsearch(@(x) sum(abs(d-G*x)),m2); % L1-norm, less emphasis on outliers
        K_ps_uw(i)=m2(2);
        if i/10000==round(i/10000)
            fprintf('%d of %d pixels processed\n',i,ps.n_ps)
        end
    end
end

ph_scla=repmat(K_ps_uw,1,size(bperp_mat,2)).*bperp_mat;

if use_small_baselines==0
    unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
    if coest_mean_vel==0
        C_ps_uw=mean(uw.ph_uw(:,unwrap_ifg_index)-ph_scla(:,unwrap_ifg_index),2);
    else
        G=[ones(length(unwrap_ifg_index),1),ps.day(unwrap_ifg_index)-ps.day(ps.master_ix)];
        m=lscov(G,[uw.ph_uw(:,unwrap_ifg_index)-ph_scla(:,unwrap_ifg_index)]',ifg_vcm(unwrap_ifg_index,unwrap_ifg_index));
        C_ps_uw=m(1,:)';
    end
else
    C_ps_uw=zeros(ps.n_ps,1);
end

oldscla=dir([sclaname,'.mat']);
if ~isempty(oldscla)
    if isfield(oldscla,'datenum')
        olddatenum=oldscla.datenum;
    else
        olddatenum=datenum(oldscla.date);
    end
    movefile([sclaname,'.mat'],['tmp_',sclaname(3:end),datestr(olddatenum,'_yyyymmdd_HHMMSS'),'.mat']);
end

stamps_save(sclaname,ph_scla,K_ps_uw,C_ps_uw,ph_ramp,ifg_vcm)

logit(1);
