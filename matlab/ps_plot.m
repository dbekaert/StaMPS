function []=ps_plot(value_type,varargin)
%function []=ps_plot(value_type,plot_flag,lims,ref_ifg,ifg_list,n_x,...
%    cbar_flag,textsize,textcolor,lon_rg,lat_rg)
% PS_PLOT plot ps values for selected ifgs
%    PS_PLOT(VALUE_TYPE,BACKGROUND,PHASE_LIMS,REF_IFG,IFG_LIST,N_X,...
%            CBAR_FLAG,TEXTSIZE,TEXTCOLOR,LON_RG,LAT_RG) 
%
%    In the case of phase, +ve values imply displacement away from the satellite
%       when the master is earlier than the slave.
%    In the case of velocities, +ve values are towards the satellite.
%
%    valid VALUE_TYPE's are:
%    'w' for wrapped phase
%    'w-d' for wrapped phase minus smoothed dem error
%    'w-o' for wrapped phase minus orbital ramps
%    'w-dm' for wrapped phase minus dem error and master AOE
%    'w-do' for wrapped phase minus dem error and orbital ramps
%    'w-dmo' for wrapped phase minus dem error, master AOE and orbital ramps
%    'p' for spatially filtered wrapped phase 
%    'u' for unwrapped phase
%    'u-d' for unwrapped phase minus dem error  
%    'u-m' for unwrapped phase minus and master AOE  
%    'u-o' for unwrapped phase minus orbital ramps 
%    'u-dm' for unwrapped phase minus dem error and master AOE  
%    'u-do' for unwrapped phase minus dem error and orbital ramps
%    'u-dmo' for unwrapped phase minus dem error, master AOE and orbital ramps
%    'u-dms' for unwrapped phase minus dem error and all AOE 
%    'u-dmos' for unwrapped phase minus dem error, all AOE and orbital ramps
%    'd' for spatially correlated DEM error (rad/m)
%    'm' for AOE phase due to master
%    'o' for orbital ramps 
%    's' for atmosphere and orbit error (AOE) phase due to slave
%    'v' mean LOS velocity (MLV) in mm/yr
%    'v-d' 
%    'v-o' 
%    'v-do' 
%    'vs' standard deviation of MLV (mm/yr)
%    'vs-d'
%    'vs-o' 
%    'vs-do' 
%    'vdrop' MLV calculated from all but current ifg (mm/yr)
%    'vdrop-d' 
%    'vdrop-o' 
%    'vdrop-do' 
%    'usb' for unwrapped phase of small baseline ifgs
%    'usb-d' 
%    'usb-o' 
%    'usb-do' 
%    'rsb' residual between unwrapped phase of sb ifgs and inverted
% 
%    When the wrapped interferograms are small baseline, 'v' and 'd' plots
%    are calulated from the unwraped small baseline interferograms by 
%    default. To force use of the single master interferograms, capitalise
%    e.g. ps_plot('V-D')
%    
%    BACKGROUND = -1 outputs the data to a .mat file instead of plotting
%                 0, black background, lon/lat axes 
%                 1, white background, lon/lat axes (default)
%                 2, shaded relief topo, lon/lat axes
%                 3, 3D topo, lon/lat axes
%                 4, mean amplitude image
%                 5, mean amplitude image, brightness showing through PS
%                 6, white background, xy axis (rotated lon/lat)
%                           
%    PHASE_LIMS = 1x2 vector with colormap limits (or 0 for default) 
%                 defaults to the range of the plotted phase 
%
%    REF_IFG = number of interferogram to reference to - defaults to 0 (master)
%              -1 for incremental referencing
%
%    IFG_LIST = list of interferograms to plot - defaults to [] (all)
%
%    N_X = maximum number of images to plot per row 
%          defaults to 0 (find optimum based on image size)
%
%    CBAR_FLAG = colorbar flag - defaults to 0 (plot on master, if plotted)
%                 1 = don't plot a colorbar
%                 2 = plot a colorbar underneath
%
%    TEXTSIZE = size of date text in points - defaults to 0 (best)
%               +ve size plots a top (default), -ve size plots at bottom
%
%    TEXTCOLOR = 1x3 color vector - default white or black depending on BACKGROUND
%
%    LON_RG = longitude range - defaults to [] (whole image)
%
%    LAT_RG = latitude range - defaults to [] (whole image)
%
%    'ts'   = produce time series plot on user click over velocity plots.
%             the position of this switch is not important.
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   09/2006 AH: 'v' option added
%   01/2007 AH: several new options added and 'u' options changed to print
%               only unwrapped IFGs by default
%   01/2007 AH: Added 64-bit machine compatibility
%   01/2009 AH: BACKGROUND = -1 option added
%   06/2009 AH: Orbital ramps added
%   09/2009 AH: Sign for velocity plots flipped
%   02/2010 AH: Give warning only if unable to save velocities
%   02/2010 AH: Replace unwrap_ifg_index with drop_ifg_index
%   03/2010 AH: Add var/cov to 'vsb' inversion and plot instead of
%   03/2010 AH: Plot 'vsb' instead of 'v' by default for small baselines
%   03/2010 AH: Rationalise velocity plotting
%   06/2010 AH: Correct bug for 'w' option for all ifgs 
%   06/2010 AH: Correct bug for 'v' option for subset small baseline ifgs
%   07/2010 DB: Correct bug for N_X plotting option
%   08/2010 DB: Option to plot phase data and edit frontsize colorbar
%   08/2010 DB: Adding figure names to plotted results
%   11/2010 AH: Fix bug for w-dmo
%   11/2010 AH: Add u-dmos
%   11/2010 MMC & MA: plotting time series using 'ts' option 
%   11/2010 DB: Changing position of colorbar for ts plot
%   ======================================================================

stdargin = nargin ;

parseplotprm  % check if 'ts' is specified

if stdargin<1
    help ps_plot
    error('not enough input args')
end
    
if stdargin < 2
    plot_flag=1;
end

if stdargin < 3
    lims=[];
end

if stdargin < 4
    ref_ifg=0;
end

if stdargin < 5
    ifg_list=[];
end

if stdargin < 6
    n_x=0;
end

if stdargin < 7
    cbar_flag=0;
end

if stdargin < 8
    textsize=0;
end

if stdargin < 9 | isempty(textcolor)
    if plot_flag==1 | plot_flag==6
        textcolor=[0 0 0.004];
    else
        textcolor=[1 1 0.996];
    end
end

if lims==0
    lims=[];
end

if stdargin < 10
    lon_rg=[];
end

if stdargin < 11
    lat_rg=[];
end

n_y=0;
units='rad';

load psver
psname=['./ps',num2str(psver)];
pmname=['./pm',num2str(psver)];
rcname=['./rc',num2str(psver)];
phuwname=['./phuw',num2str(psver)];
phuwsbname=['./phuw_sb',num2str(psver)];
phuwsbresname=['./phuw_sb_res',num2str(psver)];
scnname=['./scn',num2str(psver)];
ifgstdname=['./ifgstd',num2str(psver)];
apsname=['./aps',num2str(psver)];
apssbname=['./aps_sb',num2str(psver)];
sclaname=['./scla',num2str(psver)];
sclasbname=['./scla_sb',num2str(psver)];
sclasmoothname=['./scla_smooth',num2str(psver)];
sclasbsmoothname=['./scla_smooth_sb',num2str(psver)];
meanvname=['./mv',num2str(psver)];

ps=load(psname);
day=ps.day;
master_day=ps.master_day;
xy=ps.xy;
lonlat=ps.lonlat;
n_ps=ps.n_ps;
n_ifg=ps.n_ifg;

master_ix=sum(day<master_day)+1;

ref_ps=0;    

drop_ifg_index=getparm('drop_ifg_index');
small_baseline_flag=getparm('small_baseline_flag');

if strcmpi(small_baseline_flag,'y')
    unwrap_ifg_index_sb=setdiff([1:ps.n_ifg],drop_ifg_index);
    if value_type(1)~='w' & value_type(1)~='p'
      warning('off','MATLAB:load:variableNotFound');
      phuw=load(phuwname,'unwrap_ifg_index_sm');
      warning('on','MATLAB:load:variableNotFound');
      if isfield(phuw,'unwrap_ifg_index_sm');
        unwrap_ifg_index=phuw.unwrap_ifg_index_sm;
      else
        unwrap_ifg_index=[1:ps.n_image];
      end
    else
        unwrap_ifg_index=unwrap_ifg_index_sb;
    end

    if length(value_type)>2 & (value_type(1:3)=='usb'|value_type(1:3)=='rsb') & isempty(ifg_list)
        ifg_list=unwrap_ifg_index_sb;
    end
else
    unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);
end
if (value_type(1)=='u' | value_type(1)=='a' | value_type(1)=='w') & isempty(ifg_list)
    ifg_list=unwrap_ifg_index;
end



if ischar(value_type)==1
group_type=value_type;
if length(group_type)>1 & strcmpi(group_type(1:2),'vs')
    group_type='vs'; 
    if strcmpi(small_baseline_flag,'y') & strcmp(value_type(1:2),'vs')
        use_small_baselines=1;
        fprintf('Velocity std devs calculated from small baseline interferograms\n')
    else
        use_small_baselines=0;
    end
elseif strcmpi(group_type(1),'v')
    if strcmpi(small_baseline_flag,'y') & strcmp(group_type(1),'v')
        group_type='vsb';
        fprintf('Velocities calculated from small baseline interferograms\n')
    else
        group_type='v';
    end
elseif strcmpi(group_type(1),'d')
    if strcmpi(small_baseline_flag,'y') & strcmp(group_type(1),'d')
        group_type='dsb';
        fprintf('DEM error calculated from small baseline interferograms\n')
    else
        group_type='d';
    end
end


value_type=lower(value_type);
switch(group_type)
    case {'w'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc
	fig_name = 'w';
     case {'w-d'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        ph_all=ph_all.*exp(-j*scla.ph_scla);
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
        fig_name = 'w-d';
     case {'w-o'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclaname);
        else
            scla=load(sclasbname);
        end
        ph_all=ph_all.*exp(-j*scla.ph_ramp);
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
        fig_name = 'w-o';
     case {'w-do'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        ph_all=ph_all.*exp(-j*(scla.ph_scla+scla.ph_ramp));
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
        fig_name = 'w-do';
     case {'w-dm'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        ph_all=ph_all.*exp(-j*scla.ph_scla);
        ph_all=ph_all.*repmat(exp(-j*scla.C_ps_uw),1,size(ph_all,2));
        ph_all(:,ps.master_ix)=1;
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
        fig_name = 'w-dm';
     case {'w-dmo'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        ph_all=ph_all.*exp(-j*(scla.ph_scla+scla.ph_ramp));
        ph_all=ph_all.*repmat(exp(-j*scla.C_ps_uw),1,size(ph_all,2));
        ph_all(:,ps.master_ix)=1;
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
        fig_name = 'w-dmo';
    case {'u'}
        phuw=load(phuwname);
        ph_all=phuw.ph_uw;
        clear phuw
        ref_ps=ps_setref;
        fig_name = 'u';
    case {'usb'}
        uw=load(phuwsbname);
        ph_all=uw.ph_uw;
        clear uw
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb';
    case {'dsb'}
        scla=load(sclasbname);
        ph_all=scla.K_ps_uw;
        clear scla
        ref_ps=ps_setref;
        units='rad/m';
        fig_name = 'dsb';
    case {'usb-d'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        ph_all=uw.ph_uw - scla.ph_scla;
        clear uw scla
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb-d';
    case {'usb-o'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        ph_all=uw.ph_uw - scla.ph_ramp;
        clear uw scla
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb-o';
    case {'usb-do'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        ph_all=uw.ph_uw - scla.ph_scla - scla.ph_ramp;
        clear uw scla
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb-do';
    case {'usb-a'}
        uw=load(phuwsbname);
        aps=load(apssbname);
        ph_all=uw.ph_uw - aps.ph_aps_slave;
        clear uw aps
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb-a';
    case {'usb-da'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        aps=load(apssbname);
        ph_all=uw.ph_uw - scla.ph_scla - aps.ph_aps_slave;
        clear uw scla aps
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb-da';
    case {'rsb'}
        uw=load(phuwsbname);
        res=load(phuwsbresname);
        ph_all=zeros(size(uw.ph_uw));
        ph_all(:,unwrap_ifg_index_sb)=uw.ph_uw(:,unwrap_ifg_index_sb)-res.ph_res(:,unwrap_ifg_index_sb);
        clear uw
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'rsb';
    case {'asb'}
        aps=load(apssbname);
        ph_all=aps.ph_aps_slave;
        clear aps
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'asb';
    case {'u-dms'}
        uw=load(phuwname);
        scn=load(scnname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scn.ph_scn_slave - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla;
        clear uw scn scla
        ph_all(:,ps.master_ix)=0;
        ref_ps=ps_setref;
        fig_name = 'u-dms';
    case {'u-dmas'}
        uw=load(phuwname);
        scn=load(scnname);
        scla=load(sclaname);
        aps=load(apsname);
        ph_all=uw.ph_uw - scn.ph_scn_slave - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla - aps.ph_aps_slave;
        clear uw scn scla
        ref_ps=ps_setref;
        fig_name = 'u-dmas';
    case {'u-dm'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw; 
        %ph_all(:,unwrap_ifg_index)=ph_all(:,unwrap_ifg_index) - repmat(scla.C_ps_uw,1,length(unwrap_ifg_index));
        %ph_all(:,unwrap_ifg_index) = ph_all(:,unwrap_ifg_index) - scla.ph_scla(:,unwrap_ifg_index);
        ph_all=uw.ph_uw - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla;
        ph_all(:,master_ix)=0;
        clear uw scla
        ref_ps=ps_setref;
        fig_name = 'u-dm';
    case {'u-dmo'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla - scla.ph_ramp;
        ph_all(:,master_ix)=0;
        clear uw scla
        ref_ps=ps_setref;
        fig_name = 'u-dmo';
    case {'u-dma'}
        uw=load(phuwname);
        scla=load(sclaname);
        aps=load(apsname);
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla - aps.ph_aps_slave;
        ph_all(:,master_ix)=0;
        clear uw scla aps
        ref_ps=ps_setref;
        fig_name = 'u-dma';
    case {'u-dmos'}
        uw=load(phuwname);
        scn=load(scnname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scn.ph_scn_slave - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla - scla.ph_ramp;
        clear uw scn scla
        ph_all(:,ps.master_ix)=0;
        ref_ps=ps_setref;
        fig_name = 'u-dms';
    case {'u-a'}
        uw=load(phuwname);
        aps=load(apsname);
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - aps.ph_aps_slave;
        ph_all(:,master_ix)=0;
        clear uw aps
        ref_ps=ps_setref;
        fig_name = 'u-a';
    case {'u-d'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scla.ph_scla;
        clear uw scla
        ref_ps=ps_setref;
        fig_name = 'u-d';
    case {'u-o'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scla.ph_ramp;
        clear uw scla
        ref_ps=ps_setref;
        fig_name = 'u-o';
    case {'u-do'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scla.ph_scla - scla.ph_ramp;
        clear uw scla
        ref_ps=ps_setref;
        fig_name = 'u-do';
    case {'u-m'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw;
        ph_all=uw.ph_uw - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2));
        ph_all(:,master_ix)=0;
        clear uw scla
        ref_ps=ps_setref;
        fig_name = 'u-m';
    case {'a'}
        aps=load(apsname);
        %ph_all=exp(j*ph_scn);
        ph_all=aps.ph_aps_slave;
        clear aps
        ref_ps=ps_setref;
        fig_name = 'a';
    case {'s'}
        scn=load(scnname);
        %ph_all=exp(j*ph_scn);
        ph_all=scn.ph_scn_slave;
        clear scn
        ref_ps=ps_setref;
        fig_name = 's';
    case {'m'}
        %scn=load(scnname);
        %ph_all=scn.ph_scn_master;
        %clear scn
        scla=load(sclaname);
        ph_all=scla.C_ps_uw;
        clear scla
        ref_ps=ps_setref;
        fig_name = 'm';
    case {'d'}
        scla=load(sclaname);
        ph_all=scla.K_ps_uw;
        clear scla
        ref_ps=ps_setref;
        units='rad/m';
        fig_name = 'd';
    case {'o'}
        scla=load(sclaname);
        ph_all=scla.ph_ramp;
        clear scla
        ref_ps=ps_setref;
        units='rad/m';
        fig_name = 'o';
    case {'v'}
        uw=load(phuwname);
        ph_uw=uw.ph_uw;
        clear uw
        switch(value_type)
        case {'v'}
            fig_name = 'v';
        case {'v-d'}
            scla=load(sclaname);
            ph_uw=ph_uw - scla.ph_scla;
            clear scla
            fig_name = 'v-d';
        case {'v-o'}
            scla=load(sclaname);
            ph_uw=ph_uw - scla.ph_ramp;
            clear scla
            fig_name = 'v-o';
        case {'v-do'}
            scla=load(sclaname);
            ph_uw=ph_uw - scla.ph_ramp - scla.ph_scla;
            clear scla
            fig_name = 'v-do';
        case {'v-a'}
            aps=load(apsname);
            ph_uw=ph_uw - aps.ph_aps_slave;
            clear scla aps
            fig_name = 'v-a';
        case {'v-da'}
            scla=load(sclaname);
            aps=load(apsname);
            ph_uw=ph_uw - scla.ph_scla- aps.ph_aps_slave;
            clear scla aps
            fig_name = 'v-da';
        case {'v-ds'}
            scla=load(sclaname);
            scn=load(scnname);
            ph_uw=ph_uw - scla.ph_scla - scn.ph_scn_slave;
            clear scla scn
            fig_name = 'v-ds';
        case {'v-das'}
            scla=load(sclaname);
            aps=load(apsname);
            scn=load(scnname);
            ph_uw=ph_uw - scla.ph_scla - aps.ph_aps_slave - scn.ph_scn_slave;
            clear scla scn
            fig_name = 'v-das';
        case {'vdrop-d'}
            scla=load(sclaname);
            ph_uw=ph_uw - scla.ph_scla;
            clear scla
            fig_name = 'vdrop-d';
        otherwise
            error('unknown value type')
        end
        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        
        if unwrap_ifg_index(1)~=ps.master_ix & unwrap_ifg_index(end)~=ps.master_ix
            unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix); % need to include it if not ifgs either side of master
        end
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:),1),n_ps,1);
        % Each ifg has master APS - slave APS, including master 
        % (where slave APS = master APS) so OK to include master in inversion
        if strcmpi(small_baseline_flag,'y')
            phuwres=load(phuwsbresname,'sm_cov');
            if isfield(phuwres,'sm_cov');
                sm_cov=phuwres.sm_cov(unwrap_ifg_index,unwrap_ifg_index);
            else
                sm_cov=eye(length(unwrap_ifg_index));
            end
        else
            if ~exist([ifgstdname,'.mat',],'file')
%                ps_calc_ifg_std;
                sm_cov=eye(length(unwrap_ifg_index));
            else
                ifgstd=load(ifgstdname);
              if isfield(ifgstd,'ifg_std');
                ifgvar=(ifgstd.ifg_std*pi/181).^2;
                sm_cov=diag(ifgvar(unwrap_ifg_index));
              else
                sm_cov=eye(length(unwrap_ifg_index));
              end
            end
        end

        G=[ones(size(day)),day-master_day]; 
        lambda=getparm('lambda');

        if length(value_type)>4 & strcmpi(value_type(1:5),'vdrop') 
            ph_all=zeros(size(ph_uw));
            n=size(ph_uw,2);
            for i=1:n
                m=lscov(G([1:i-1,i+1:end],:),double(ph_uw(:,[1:i-1,i+1:n])',sm_cov));
                ph_all(:,i)=-m(2,:)'*365.25/4/pi*lambda*1000; 
            end
        else 
            m=lscov(G,double(ph_uw'),sm_cov);
            ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        end

        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        textsize=0;
        units='mm/yr';
        
        % TS PLOT preparation
        if ts_flag==1
           ts_flaghelper; % script to save neccessary matrices for TS plot
        end
        
   case {'vsb'}
        phuw=load(phuwsbname);
        ref_ps=ps_setref;
        ph_uw=phuw.ph_uw;
        clear phuw
        switch(value_type)
        case('v')
            fig_name = 'v';
        case('v-d')
            scla=load(sclasbname);
            ph_uw=ph_uw - scla.ph_scla;
            clear scla
            fig_name = 'v-d';
        case('v-o')
            scla=load(sclasbname);
            ph_uw=ph_uw - scla.ph_ramp;
            clear scla
            fig_name = 'v-o';
        case('v-do')
            scla=load(sclasbname);
            ph_uw=ph_uw - scla.ph_scla - scla.ph_ramp;
            clear scla
            fig_name = 'v-do';
        case {'vdrop-d'}
            scla=load(sclasbname);
            ph_uw=ph_uw - scla.ph_scla;
            clear scla
            fig_name = 'vdrop-d';
        otherwise
            error('unknown value type')
        end
        
        if ~isempty(ifg_list)
            unwrap_ifg_index_sb=intersect(unwrap_ifg_index_sb,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index_sb);
        phuwres=load(phuwsbresname,'sb_cov');
        if isfield(phuwres,'sb_cov');
            sb_cov=phuwres.sb_cov(unwrap_ifg_index_sb,unwrap_ifg_index_sb);
        else
            sb_cov=eye(length(unwrap_ifg_index_sb));
        end
        ifgday_ix=ps.ifgday_ix(unwrap_ifg_index_sb,:);
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:),1),n_ps,1);
        G=[ones(size(ifgday_ix(:,1))),day(ifgday_ix(:,2))-day(ifgday_ix(:,1))];
        lambda=getparm('lambda');

        if length(value_type)>4 & strcmpi(value_type(1:5),'vdrop') 
            ph_all=zeros(size(ph_uw));
            n=size(ph_uw,2);
            for i=1:n
                m=lscov(G([1:i-1,i+1:end],:),double(ph_uw(:,[1:i-1,i+1:n])'),sb_cov);
                ph_all(:,i)=-m(2,:)'*365.25/4/pi*lambda*1000; 
            end
        else
            m=lscov(G,double(ph_uw'),sb_cov);
            ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; 
        end
        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        textsize=0;
        units='mm/yr';
    case {'p'}
        pm=load(pmname);
        ph_all=pm.ph_patch./abs(pm.ph_patch);
        if n_ifg~=size(ph_all,2)
            ph_all=[ph_all(:,1:ps.master_ix-1),zeros(ps.n_ps,1),ph_all(:,ps.master_ix:end)];
        end
        clear pm
        if ref_ifg~=0
            ph_all=ph_all.*repmat(conj(ph_all(:,ref_ifg)),1,n_ifg);
            ph_all(:,ref_ifg)=1; % may not be so because of rounding errors
        end
        fig_name = 'p';
    case {'wf'}
        uw=load('uw_grid');
        gridix=zeros(size(uw.nzix));
        gridix(uw.nzix)=[1:uw.n_ps];
        ph_all=zeros(ps.n_ps,uw.n_ifg);
        for i=1:ps.n_ps
            ph_all(i,:)=uw.ph(gridix(uw.grid_ij(i,1),uw.grid_ij(i,2)),:);
        end
        clear uw
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(ph_all(:,ref_ifg),1,n_ifg));
        end
        fig_name = 'wf';
    case {'vs'}
        ps_mean_v(ifg_list,200,value_type(4:end),use_small_baselines);
        ifg_list=[];
        mv=load(meanvname);
        ph_all=mv.mean_v_std;
        units='mm/yr';
        fig_name = 'vs';
    otherwise
        error('unknown value type')
    end
else
	ph_all = value_type;
	fig_name = 'data';
end


if isempty(ifg_list)
    ifg_list=1:size(ph_all,2);
end
n_ifg_plot=length(ifg_list);


xgap=0.1;
ygap=0.2;
[Y,X]=meshgrid([0.7:-1*ygap:0.1],[0.1:xgap:0.8]);

if ~isempty(lon_rg)
    ix=lonlat(:,1)>=lon_rg(1)&lonlat(:,1)<=lon_rg(2);
    lonlat=lonlat(ix,:);
end

if ~isempty(lat_rg)
    ix=lonlat(:,2)>=lat_rg(1)&lonlat(:,2)<=lat_rg(2);
    lonlat=lonlat(ix,:);
end

max_xy=llh2local([max(lonlat),0]',[min(lonlat),0]);
fig_ar=4/3; % aspect ratio of figure window
useratio=1; % max fraction of figure window to use
n_i=max_xy(2)*1000;
n_j=max_xy(1)*1000;
ar=max_xy(1)/max_xy(2); % aspect ratio (x/y)
if n_x==0
    n_y=ceil(sqrt((n_ifg_plot)*ar/fig_ar)); % number of plots in y direction
    n_x=ceil((n_ifg_plot)/n_y);
    fixed_fig = 1; 			% figure with fixed aspect ratio
else
    n_y=ceil((n_ifg_plot)/n_x);
    fixed_fig = 0;
end


d_x=useratio/n_x;
d_y=d_x/ar*fig_ar;
if d_y>useratio/n_y & fixed_fig==1	% TS figure exceeds fig size
    d_y=useratio/n_y; 
    d_x=d_y*ar/fig_ar;
    h_y=0.95*d_y;
    h_x=h_y*ar/fig_ar;

    fig_size=0;
elseif d_y>useratio/n_y & fixed_fig==0 
    h_y=0.95*d_y;
    h_x=h_y*ar/fig_ar;

    y_scale = d_y*n_y;
    d_y=d_y/y_scale;   
    fig_size=1;		% check to indicate fig needs to be adapted
    h_y=0.95*d_y;

else
    h_y=0.95*d_y;
    h_x=h_y*ar/fig_ar;
    fig_size=0;
end
y=1-d_y:-d_y:0;
x=1-useratio:d_x:1-d_x;


[imY,imX]=meshgrid(y,x);
if textsize==0
    textsize=round(10*4/n_x);
    if textsize>16
        textsize=16;
    elseif textsize<8
        textsize=8;
    end
end

l_t=1/9*abs(textsize)/10; % text length
h_t=1/50*abs(textsize)/10; % text height
x_t=round((h_x-l_t)/h_x/2*n_j);
y_t=round(h_t*1.2/h_y*n_i);



ph_disp=ph_all(:,ifg_list);

if isreal(ph_all)

    if ref_ifg~=0
        if ref_ifg==-1
            ph_disp=ph_disp-[ph_disp(:,1),ph_disp(:,1:end-1)];
        else
            ph_disp=ph_disp-repmat(ph_all(:,ref_ifg),1,size(ph_disp,2));
        end
    else
        ref_ifg=master_ix;
    end
    if ref_ps~=0
        ref_ph=(ph_disp(ref_ps,:));
        mean_ph=zeros(1,size(ph_disp,2));
        for i=1:size(ph_disp,2)
            mean_ph(i)=mean(ref_ph(~isnan(ref_ph(:,i)),i),1);
        end
        ph_disp=ph_disp-repmat(mean_ph,n_ps,1);
    end
    
    phsort=sort(ph_disp(~isnan(ph_disp)));
    if isempty(lims)
        maxph=phsort(round(length(phsort)*.999));
        minph=phsort(ceil(length(phsort)*.001));
       lims=[minph,maxph];
    end
else
    if ref_ifg==0
        ref_ifg=master_ix;
    elseif ref_ifg==-1
        ph_disp=ph_disp.*conj([ph_disp(:,1),ph_disp(:,1:end-1)]);
    end
    if ref_ps~=0
        ph_disp=ph_disp./abs(ph_disp);
        ref_ph=(ph_disp(ref_ps,:));
        mean_ph=zeros(1,size(ph_disp,2));
        for i=1:size(ph_disp,2)
            mean_ph(i)=sum(ref_ph(~isnan(ref_ph(:,i)),i));
        end
        ph_disp=ph_disp.*conj(repmat(mean_ph,n_ps,1));
    end
    lims=[-pi,pi];

end

if plot_flag==-1
    savename=['ps_plot_',value_type];
    save(savename,'ph_disp','ifg_list')
else
  h_fig = figure;
  set(gcf,'renderer','zbuffer','name',fig_name)

  if fig_size==1
      Position = get(gcf,'Position');
      Position(1,2)=50;
      Position(1,4)=Position(1,2)+Position(1,4)*y_scale;
      set(gcf,'Position',Position)
  end


  i_im=0;
  for i=ifg_list
    i_im=i_im+1;
    if n_ifg_plot>1
        axes('position',[imX(i_im),imY(i_im),h_x,h_y])
    end
    ps_plot_ifg(ph_disp(:,i_im),plot_flag,lims,lon_rg,lat_rg);
    %plot_phase(ph_tc(:,i)*conj(ph_tc(ref_ix,i)));
    box on
    if n_ifg_plot>1
        set(gca,'yticklabel',[])
        set(gca,'xticklabel',[])
    end
    xlim=get(gca,'xlim');
    x_t=(h_x-l_t)/2/h_x*(xlim(2)-xlim(1))+xlim(1);
    ylim=get(gca,'ylim');
    if textsize>0
        y_t=(h_y-1.2*h_t)/h_y*(ylim(2)-ylim(1))+ylim(1);
    else
        y_t=(0.5*h_t)/h_y*(ylim(2)-ylim(1))+ylim(1);
    end
    %xlabel([num2str((day(i)/365.25),3),'yr, ',num2str(round(bperp(i))),'m'])
    if textsize~=0 & size(day,1)==size(ph_all,2)
        t=text(x_t,y_t,[datestr(day(i),'dd mmm yyyy')]);
        set(t,'fontweight','bold','color',textcolor,'fontsize',abs(textsize))
    end
    
    if cbar_flag==0 & (i==ref_ifg | (isempty(intersect(ref_ifg,ifg_list)) & i==ifg_list(1))) 
        if n_ifg_plot>1
            h=colorbar('South');
	    xlim=get(h,'xlim');
	    set(h,'xlim',[xlim(2)-64,xlim(2)])

        else
            %h=colorbar('SouthOutside');
	    h = colorbar('peer',gca);
            ylim=get(h,'ylim');
            set(h,'ylim',[ylim(2)-64,ylim(2)])
        end

        if diff(lims)>1 | diff(lims)==0
            plotlims=round(lims*10)/10;
        else
            limorder=ceil(-log10(diff(lims)))+2;
            plotlims=round(lims*10^limorder)/10^limorder;
        end
	if n_ifg_plot>1
        	set(h,'xtick',[xlim(2)-64,xlim(2)],'Xticklabel',plotlims,'xcolor',textcolor,'ycolor',textcolor,'fontweight','bold','color',textcolor,'FontSize',abs(textsize))
	        h=xlabel(h,units);
        	pos=get(h,'position');
        	pos(2)=pos(2)/2.2;
        	set(h,'position',pos,'FontSize',abs(textsize));
	else
                set(h,'ytick',[ylim(2)-64,ylim(2)],'yticklabel',plotlims,'xcolor',textcolor,'ycolor',textcolor,'fontweight','bold','color',textcolor,'FontSize',abs(textsize))
	        set(get(h,'ylabel'),'String',units,'FontSize',abs(textsize))  

	end
    end

  end
end
    

fprintf('Color Range: %g to %g %s\n',lims,units)

if ts_flag == 1
  figure(h_fig);
  clear all % clean up to save memory
  
  % Place button for new TS plots
  fPosition=get(gcf,'position');
  %  pos_new = [pos(1) pos(2)+0.1 pos(3) pos(4)-0.1]
  
  % new TS plot button
  mButton=uicontrol('Style', 'pushbutton', 'Callback', 'ts_plot',...
    'String','new TS plot', 'Position', [150 30 90 20] ); % if exist don't create
  %                                      left bottom width height
  %mButtonposition=get(mButton,'Position')
 
  % Radius Text boxes
  mTextBox=uicontrol('Style', 'edit','String','radius = 2.778e-4 deg. x', 'Position',...
      [260 30 150 20] );

  mEditBox=uicontrol('Style', 'edit','String','2', 'Position',...
      [410 30 30 20],'BackgroundColor',[1 1 1] );
  ts_plot   % select a point then plot
end


    
    

