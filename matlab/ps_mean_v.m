function [fig_name_tca]=ps_mean_v(ifg_list,n_boot,subtract_switches,use_small_baselines,aps_flag)
%PS_MEAN_V Calculate mean velocities and their standard deviations
%   PS_MEAN_V(IFG_LIST,N_BOOT,SUBTRACT_SWITCHES) where IFG_LIST gives indices 
%   of interferograms to be included in the calculation,N_BOOT specifies number 
%   of bootstrap iterations and RAMP_FLAG is 1 to subtract orbital ramps
%
%   Andy Hooper, Jan 2008
%
%   ======================================================================
%   06/2009 AH: Subtract orbital ramp if present
%   06/2009 AH: Process smaller chunks to reduce memory needs
%   02/2010 AH: Revert to v3.1 version
%   02/2010 AH: Replace unwrap_ifg_index with drop_ifg_index
%   03/2010 AH: 3rd input changed to 'subtract_switches'
%   03/2010 AH: Use small baselines option added
%   03/2010 AH: var/cov added to inversion
%   05/2010 AH: Include master if there are not ifgs before and after
%   06/2010 AH: small change to ps_mean_v.m
%   03/2014 AH: -a etc added
%   04/2014 DB: fix fig_name_tca output in case no aps correction
%   05/2014 DB: Fix v-doa to v-dao 
%   05/2014 DB: For APS related corrections include deramping on the fly
%   05/2014 DB: big fix in case of a nan 
%   07/2014 EH: Fix for nanmean in case of single element compared to vector
%   03/2015 DB: Do all deramping on the fly, make consistent for TRAIN release
%   05/2017 DB: use stamps save to save larger variables than 2GB.
%   06/2017 DB: Include the option to choose between chol decompostion for
%               bootstrapping (not to invert covariance often), or conventional
%               bootstrapping requiring to invert covariance matrix for each iteration.
%   ======================================================================


fprintf('Calculating standard deviation of mean velocity...\n')


chol_flag = 'n';

if nargin<1
    ifg_list=[];
end
if nargin<2
    n_boot=100;
end

if nargin<3 | isempty(subtract_switches)
   subtract_switches='';
end

if nargin<4
    use_small_baselines=0;
end

if nargin<5
    aps_flag=1;
end

load psver
psname=['./ps',num2str(psver)];
phuwname=['./phuw',num2str(psver)];
phuwsbresname=['./phuw_sb_res',num2str(psver)];
ifgstdname=['./ifgstd',num2str(psver)];
sclaname=['./scla',num2str(psver)];
apsname=['./tca',num2str(psver)];
tidename=['./tide',num2str(psver)];

mvname=['mv',num2str(psver)];

ps=load(psname);

drop_ifg_index=getparm('drop_ifg_index');
fig_name_tca = '';
if strcmpi(getparm('small_baseline_flag'),'y')
    if use_small_baselines==0
        phuw=load(phuwname,'unwrap_ifg_index_sm');
        if isfield(phuw,'unwrap_ifg_index_sm');
            unwrap_ifg_index=phuw.unwrap_ifg_index_sm;
        else
            unwrap_ifg_index=[1:ps.n_image];
        end
        phuwres=load(phuwsbresname,'sm_cov');
        if isfield(phuwres,'sm_cov');
            ifg_cov=phuwres.sm_cov;
        else
            ifg_cov=eye(ps.n_image);
        end
    else
        unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);
        phuwname=['./phuw_sb',num2str(psver)];
        sclaname=['./scla_sb',num2str(psver)];
        apsname=['./tca_sb',num2str(psver)];
        tidename=['./tide_sb',num2str(psver)];

        phuwres=load(phuwsbresname,'sb_cov');
        if isfield(phuwres,'sb_cov');
            ifg_cov=phuwres.sb_cov;
        else
            ifg_cov=eye(ps.n_ifg);
        end
    end
else
    use_small_baselines=0;
    unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);
    if ~exist([ifgstdname,'.mat'],'file')
        ifg_cov=eye(ps.n_ifg);
    else ifgstd=load(ifgstdname);
      if isfield(ifgstd,'ifg_std');
        ifgvar=(ifgstd.ifg_std*pi/181).^2;
        ifg_cov=diag(ifgvar);
      else
        ifg_cov=eye(ps.n_ifg);
      end
    end
end



uw=load(phuwname);
ph_uw=uw.ph_uw;
clear uw


switch(subtract_switches)
case('')
case('d')
    scla=load(sclaname);
    ph_uw=ph_uw - scla.ph_scla;
    clear scla
case('o')
    [ph_uw] = ps_deramp(ps,ph_uw);
case('a')
    aps=load(apsname);
    [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
    ph_uw=ph_uw - aps_corr;
    clear aps aps_corr
case('ao')
    aps=load(apsname);
    [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
    ph_uw=ph_uw - aps_corr;
    [ph_uw] = ps_deramp(ps,ph_uw);
    clear aps aps_corr
case('do')
    scla=load(sclaname);
    ph_uw=ph_uw - scla.ph_scla;
    [ph_uw] = ps_deramp(ps,ph_uw);
case('da')
    scla=load(sclaname);
    ph_uw=ph_uw - scla.ph_scla;
    clear scla
    aps=load(apsname);
    [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
    ph_uw=ph_uw - aps_corr;
    clear aps aps_corr
case('dat')
    tide=load(tidename);
    scla=load(sclaname);
    ph_uw=ph_uw - scla.ph_scla;
    clear scla
    aps=load(apsname);
    [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
    ph_uw=ph_uw - aps_corr-tide.ph_tide;
    clear aps aps_corr
case('dao')
    scla=load(sclaname);
    ph_uw=ph_uw - scla.ph_scla ;
    [ph_uw] = ps_deramp(ps,ph_uw);
    aps=load(apsname);
    [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
    ph_uw=ph_uw - aps_corr;
    clear aps aps_corr
otherwise
        error('unknown subtract flags')
end

if use_small_baselines==0 
    %AH2CHECK %%%%%%% & unwrap_ifg_index(1)~=ps.master_ix & unwrap_ifg_index(end)~=ps.master_ix
    unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
end

ph_all=zeros(ps.n_ps,1);
ref_ps=ps_setref;
if ~isempty(ifg_list)
    unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
    ifg_list=[];
end
ph_uw=ph_uw(:,unwrap_ifg_index);
ifg_cov=ifg_cov(unwrap_ifg_index,unwrap_ifg_index);
if rank(ifg_cov)<size(ifg_cov,1)
    error('Interferograms included for which there is no unwrapped solution')
end

if use_small_baselines==0
    day=ps.day(unwrap_ifg_index)-ps.master_day;
else 
    day=ps.day(ps.ifgday_ix(unwrap_ifg_index,2))-ps.day(ps.ifgday_ix(unwrap_ifg_index,1));
end
N=length(day);
ph_uw=ph_uw-repmat(nanmean(ph_uw(ref_ps,:),1),ps.n_ps,1);
lambda=getparm('lambda');
ph_uw=double(ph_uw/4/pi*lambda*1000)'; 
G=[ones(N,1),double(day)/365.25] ;

% computing mean velocity
mean_v=lscov(G,ph_uw,ifg_cov)';
mean_v=mean_v(:,2);


% perform the bootstrapping
n=10000; % process n PS at a time
i=1;
mean_v_std=NaN(ps.n_ps,1,'single');
rand_init=sum(100*clock);
n_ifg = size(ph_uw,1);


if strcmpi(chol_flag,'y')
    % cholesky decomposition of the set of linear equations including
    % covariance. This allows to use the decomposition values in a uniform
    % LSQ inversion. i.e. the covariance matrix does not need to be
    % inverted each time.
    ph_uw=chol(inv(ifg_cov))*ph_uw;
    G=chol(inv(ifg_cov))*G;

    % looping over the PS segments
    while i<ps.n_ps
        if i+n< ps.n_ps
           i_end=i+n-1;
        else
           i_end=ps.n_ps;
        end
        ph_bit=ph_uw(:,i:i_end);

        % in-case of nans and only a subset being processed
        ix_non_nan = sum(isnan(ph_bit),1)==0;
        
        % fix the random generator, such its the same for each segment
        rand('twister',rand_init);

        % skip the segment in case of nan's
        if sum(ix_non_nan)>0
            if sum(ix_non_nan)== length(ix_non_nan)
                [mean_v_dist,boot_ix] = bootstrp(n_boot, @(x) single(lscov(G(x,:),ph_bit(x,:))), [1:N]);
                 temp = std(mean_v_dist(:,2:2:end))';
            else
                [mean_v_dist_temp,boot_ix] = bootstrp(n_boot, @(x) single(lscov(G(x,:),ph_bit(x,ix_non_nan))), [1:N]);
                temp = NaN([size(ph_bit,2) 1]); 
                temp(ix_non_nan) = std(mean_v_dist_temp(:,2:2:end))';
                clear mean_v_dist_temp
            end
            mean_v_std(i:i_end) =temp;       
        end
        i=i+n;
        fprintf('%d PS processed\n',i-1)
    end
    
    
else
    % conventional bootstrapping where the coveriance matrix is inverted on
    % each set PS segments.
    while i<ps.n_ps
        if i+n< ps.n_ps
           i_end=i+n-1;
        else
           i_end=ps.n_ps;
        end
        ph_bit=ph_uw(:,i:i_end);

        % in-case of nans and only a subset being processed
        ix_non_nan = sum(isnan(ph_bit),1)==0;
        
        % fix the random generator, such its the same for each segment
        rand('twister',rand_init)

        % skip the segment in case of nan's
        if sum(ix_non_nan)>0
            if sum(ix_non_nan)== length(ix_non_nan)
                [mean_v_dist,boot_ix] = bootstrp(n_boot, @(x) single(lscov(G(x,:),ph_bit(x,:),ifg_cov)), [1:N]);
                temp = std(mean_v_dist(:,2:2:end))';
            else
               fprintf('test\n')
               [mean_v_dist_temp,boot_ix] = bootstrp(n_boot, @(x) single(lscov(G(x,:),ph_bit(x,ix_non_nan),ifg_cov(ix_non_nan,ix_non_nan))), [1:N]);
               temp = NaN([size(ph_bit,2) 1]); 
               temp(ix_non_nan) = std(mean_v_dist_temp(:,2:2:end))';
               clear mean_v_dist_temp
            end
            mean_v_std(i:i_end) =temp;       
        end
        i=i+n;
        fprintf('%d PS processed\n',i-1)
    end
end

stamps_save(mvname,n_boot,subtract_switches,mean_v,mean_v_std)

