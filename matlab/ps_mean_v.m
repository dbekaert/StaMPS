function []=ps_mean_v(ifg_list,n_boot,ramp_flag)
%PS_MEAN_V Calculate mean velocities and their standard deviations
%   PS_MEAN_V(IFG_LIST,N_BOOT,RAMP_FLAG) where IFG_LIST gives indices of 
%   interferograms to be included in the calculation,N_BOOT specifies number 
%   of bootstrap iterations and RAMP_FLAG is 1 to subtract orbital ramps
%
%   Andy Hooper, Jan 2008
%
%   ======================================================================
%   06/2009 AH: Subtract orbital ramp if present
%   06/2009 AH: Process smaller chunks to reduce memory needs
%   02/2010 AH: Revert to v3.1 version
%   ======================================================================


fprintf('Calculating standard deviation of mean velocity...\n')

if nargin<1
    ifg_list=[];
end

if nargin<2
    n_boot=100;
end

if nargin<3
    ramp_flag=0;
end

load psver
psname=['ps',num2str(psver)];
phuwname=['phuw',num2str(psver)];
sclaname=['scla',num2str(psver)];
mvname=['mv',num2str(psver)];

ps=load(psname);

unwrap_ifg_index=getparm('unwrap_ifg_index');

if strcmpi(getparm('small_baseline_flag'),'y')
    unwrap_ifg_index=[1:ps.n_image];
else
    if strcmp(unwrap_ifg_index,'all') 
        unwrap_ifg_index=[1:ps.n_ifg];
    end
end


uw=load(phuwname);
scla=load(sclaname);
ph_uw=uw.ph_uw - scla.ph_scla;
if ramp_flag ~=0
    if isfield(scla,'ph_ramp') & size(scla.ph_ramp,1)==ps.n_ps
        ph_uw=uw.ph_uw - scla.ph_ramp;
    else
        error(['ph_ramp not present or wrong size in ',sclaname])
    end
end

clear uw scla
%unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
ph_all=zeros(ps.n_ps,1);
ref_ps=ps_setref;
if ~isempty(ifg_list)
    unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
    ifg_list=[];
end
ph_uw=ph_uw(:,unwrap_ifg_index);
day=ps.day(unwrap_ifg_index);
N=length(day);
ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),ps.n_ps,1);
lambda=getparm('lambda');
ph_uw=double(ph_uw/4/pi*lambda*1000)'; 
G=[ones(N,1),double(day-ps.master_day)/365.25] ;

mean_v=(G\ph_uw)';
mean_v=mean_v(:,2);

n=10000; % process n PS at a time
i=1;
mean_v_std=zeros(ps.n_ps,1,'single');
rand_init=sum(100*clock);

while i<ps.n_ps
    if i+n< ps.n_ps
       i_end=i+n-1;
    else
       i_end=ps.n_ps;
    end
    ph_bit=ph_uw(:,i:i_end);
    rand('twister',rand_init)
    [mean_v_dist,boot_ix] = bootstrp(n_boot, @(x) single(G(x,:)\ph_bit(x,:)), [1:N]);
    mean_v_std(i:i_end) = std(mean_v_dist(:,2:2:end))';
    i=i+n;
    fprintf('%d PS processed\n',i-1)
end

save(mvname,'n_boot','ramp_flag','mean_v','mean_v_std')

