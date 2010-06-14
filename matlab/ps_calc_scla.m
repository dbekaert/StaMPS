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
%   ================================================================
logit;
fprintf('Estimating spatially-correlated look angle error...\n')

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

if use_small_baselines~=0
    if small_baseline_flag~='y'
        error('   Use small baselines requested but there are none')
    end
end

if use_small_baselines==0
    recalc_index=getparm('recalc_index',1);
else
    recalc_index=getparm('sb_recalc_index',1);
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
    apsname=['./aps',num2str(psver)];
else
    phuwname=['./phuw_sb',num2str(psver)];
    sclaname=['./scla_sb',num2str(psver)];
    apsname=['./aps_sb',num2str(psver)];
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

if strcmpi(scla_deramp,'y')
    fprintf('\n   deramping ifgs...\n')
    ph_ramp=zeros(ps.n_ps,n_ifg,'single');
    G=double([ones(ps.n_ps,1),ps.xy(:,2),ps.xy(:,3)]);
    for i=1:length(unwrap_ifg_index)
        d=uw.ph_uw(:,unwrap_ifg_index(i));
        m=G\double(d(:));
        ph_this_ramp=G*m;
        uw.ph_uw(:,unwrap_ifg_index(i))=uw.ph_uw(:,unwrap_ifg_index(i))-ph_this_ramp; % subtract ramp
        ph_ramp(:,unwrap_ifg_index(i))=ph_this_ramp;
    end
else
    ph_ramp=[];
end


if ~strcmpi(recalc_index,'all')
    unwrap_ifg_index=intersect(unwrap_ifg_index,recalc_index);
end

if exist([apsname,'.mat'],'file')
    aps=load(apsname);
    uw.ph_uw=uw.ph_uw-aps.ph_aps_slave;
end

ref_ps=ps_setref;
uw.ph_uw=uw.ph_uw-repmat(mean(uw.ph_uw(ref_ps,:)),ps.n_ps,1);

if use_small_baselines==0
    if strcmpi(small_baseline_flag,'y')
        bperp_mat=zeros(ps.n_ps,ps.n_image,'single');
        G=zeros(ps.n_ifg,ps.n_image);
        for i=1:ps.n_ifg
             G(i,ps.ifgday_ix(i,1))=-1;
             G(i,ps.ifgday_ix(i,2))=1;
        end
        if isfield(uw,'unwrap_ifg_index_sm')
            unwrap_ifg_index=setdiff(uw.unwrap_ifg_index_sm,ps.master_ix)
        else
            unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix)
        end
        G=G(:,unwrap_ifg_index);
        bperp_some=[G\double(bp.bperp_mat')]';
        bperp_mat(:,unwrap_ifg_index)=bperp_some;
        clear bperp_some
    else

        bperp_mat=[bp.bperp_mat(:,1:ps.master_ix-1),zeros(ps.n_ps,1,'single'),bp.bperp_mat(:,ps.master_ix:end)];
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
    end
    day=ps.day-ps.master_day;

else
    bperp_mat=bp.bperp_mat;
    day=ps.ifgday(:,2)-ps.ifgday(:,1);
end

bperp=mean(bperp_mat);
fprintf('\n%d ifgs used in estimation:\n',length(unwrap_ifg_index))
for i=1:length(unwrap_ifg_index)
    if use_small_baselines~=0
        fprintf('   %s to %s %5d days %5d m\n',datestr(ps.ifgday(unwrap_ifg_index(i),1)),datestr(ps.ifgday(unwrap_ifg_index(i),2)),day(unwrap_ifg_index(i)),round(bperp(unwrap_ifg_index(i))))
    else
        fprintf('   %s %5d days %5d m\n',datestr(ps.day(unwrap_ifg_index(i))),day(unwrap_ifg_index(i)),round(bperp(unwrap_ifg_index(i))))
    end
end
fprintf('\n')

K_ps_uw=zeros(ps.n_ps,1);
C_ps_uw=zeros(ps.n_ps,1);

if coest_mean_vel==0 | length(unwrap_ifg_index)<4
    G=[ones(length(unwrap_ifg_index),1),double(bperp_mat(i,unwrap_ifg_index)')];
    m0=[0;0];
else
    G=[ones(length(unwrap_ifg_index),1),double(bperp_mat(i,unwrap_ifg_index)'),double(day(unwrap_ifg_index))];
    m0=[0;0;0];
end

ifg_vcm=eye(size(uw.ph_uw,2));
    
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
        ifg_vcm=diag((ifgstd.ifg_std*pi/180).^2);
        clear ifgstd
    end
end

if strcmpi(scla_method,'L1')
    for i=1:ps.n_ps
        m=fminsearch(@(x) sum(abs(d-G*x)),m); % L1-norm, less emphasis on outliers
        K_ps_uw(i)=m(2);
        if use_small_baselines==0
            C_ps_uw(i)=m(1);
        end
        if i/100000==round(i/100000)
            fprintf('%d of %d pixels processed\n',i,ps.n_ps)
        end
    end
else
    d=double(uw.ph_uw(:,unwrap_ifg_index)');
    m=lscov(G,d,ifg_vcm(unwrap_ifg_index,unwrap_ifg_index)); % L2-norm
    K_ps_uw=m(2,:)';
    if use_small_baselines==0
        C_ps_uw=m(1,:)';
    end
end

ph_scla=repmat(K_ps_uw,1,size(bperp_mat,2)).*bperp_mat;

oldscla=dir([sclaname,'.mat']);
if ~isempty(oldscla)
    if isfield(oldscla,'datenum')
        olddatenum=oldscla.datenum;
    else
        olddatenum=datenum(oldscla.date);
    end
    movefile([sclaname,'.mat'],['tmp_',sclaname(3:end),datestr(olddatenum,'_yyyymmdd_HHMMSS'),'.mat']);
end

save(sclaname,'ph_scla','K_ps_uw','C_ps_uw','ph_ramp','ifg_vcm')
    
logit(1);
