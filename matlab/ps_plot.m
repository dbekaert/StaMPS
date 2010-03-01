function []=ps_plot(value_type,plot_flag,lims,ref_ifg,ifg_list,n_x,cbar_flag,textsize,textcolor,lon_rg,lat_rg)
% PS_PLOT plot ps values for selected ifgs
%    PS_PLOT(VALUE_TYPE,BACKGROUND,PHASE_LIMS,REF_IFG,IFG_LIST,N_X,CBAR_FLAG,TEXTSIZE,TEXTCOLOR,LON_RG,LAT_RG) 
%
%    In the case of phase, +ve values imply displacement away from the satellite
%       when the master is earlier than the slave.
%    In the case of velocities, +ve values are towards the satellite.
%
%    valid VALUE_TYPE's are:
%    'w' for wrapped phase
%    'p' for spatially filtered wrapped phase 
%    'u' for unwrapped phase
%    'd' for spatially correlated DEM error (rad/m)
%    'm' for AOE phase due to master
%    'o' for orbital ramps 
%    's' for atmosphere and orbit error (AOE) phase due to slave
%    'w-d' for wrapped phase minus smoothed dem error
%    'w-o' for wrapped phase minus orbital ramps
%    'w-dm' for wrapped phase minus dem error and master AOE
%    'w-do' for wrapped phase minus dem error and orbital ramps
%    'w-dmo' for wrapped phase minus dem error, master AOE and orbital ramps
%    'u-d' for unwrapped phase minus dem error
%    'u-m' for unwrapped phase minus master AOE
%    'u-o' for unwrapped phase minus orbital ramps
%    'u-dm' for unwrapped phase minus dem error and master AOE
%    'u-do' for unwrapped phase minus dem error and orbital ramps
%    'u-dmo' for unwrapped phase minus dem error, master AOE and orbital ramps
%    'u-dms' for unwrapped phase minus dem error and all AOE
%    'v' mean LOS velocity (MLV) in mm/yr
%    'v-d' MLV calculated after removal of dem error (mm/yr)
%    'v-o' MLV calculated after removal of orbital ramps (mm/yr)
%    'v-do' MLV calculated after removal of dem error and orbital ramps (mm/yr)
%    'vs-d' standard deviation of MLV after removal of dem error (mm/yr)
%    'vs-do' standard deviation of MLV after removal of dem error and orbital ramps (mm/yr)
%    'usb' for unwrapped phase of small baseline ifgs
%    'dsb' for spatially correlated DEM error (rad/m) 
%    'usb-d' for unwrapped phase of sb ifgs minus dem error 
%    'usb-o' for unwrapped phase of sb ifgs minus orbital ramps
%    'usb-do' for unwrapped phase of sb ifgs minus dem error and orbital ramps
%    'rsb' residual between unwrapped phase of sb ifgs and inverted
%    'vsb' mean LOS velocity (mm/yr) from sb ifgs
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
%   ======================================================================



if nargin<1
    help ps_plot
    error('not enough input args')
end
    
if nargin < 2
    plot_flag=1;
end

if nargin < 3
    lims=[];
end

if nargin < 4
    ref_ifg=0;
end

if nargin < 5
    ifg_list=[];
end

if nargin < 6
    n_x=0;
end

if nargin < 7
    cbar_flag=0;
end

if nargin < 8
    textsize=0;
end

if nargin < 9 | isempty(textcolor)
    if plot_flag==1 | plot_flag==6
        textcolor=[0 0 0.004];
    else
        textcolor=[1 1 0.996];
    end
end

if lims==0
    lims=[];
end

if nargin < 10
    lon_rg=[];
end

if nargin < 11
    lat_rg=[];
end

n_y=0;
units='rad';

load psver
psname=['ps',num2str(psver)];
pmname=['pm',num2str(psver)];
rcname=['rc',num2str(psver)];
rcuwname=['rcuw',num2str(psver)];
phuwname=['phuw',num2str(psver)];
phuwsbname=['phuw_sb',num2str(psver)];
phuwsbresname=['phuw_sb_res',num2str(psver)];
scnname=['scn',num2str(psver)];
apsname=['aps',num2str(psver)];
apssbname=['aps_sb',num2str(psver)];
sclaname=['scla',num2str(psver)];
sclasbname=['scla_sb',num2str(psver)];
sclasmoothname=['scla_smooth',num2str(psver)];
sclasbsmoothname=['scla_smooth_sb',num2str(psver)];
meanvname=['mv',num2str(psver)];

ps=load(psname);
day=ps.day;
master_day=ps.master_day;
xy=ps.xy;
lonlat=ps.lonlat;
n_ps=ps.n_ps;
n_ifg=ps.n_ifg;

master_ix=sum(day<master_day)+1;

ref_ps=0;    

unwrap_ifg_index=getparm('unwrap_ifg_index');

if strcmpi(getparm('small_baseline_flag'),'y')
    unwrap_ifg_index_sb=unwrap_ifg_index;
    if strcmpi(unwrap_ifg_index,'all')
        unwrap_ifg_index_sb=[1:ps.n_ifg];
    end
    
    if value_type(1)~='w' & value_type(1)~='p'
      warning('off','MATLAB:load:variableNotFound');
      phuw=load(phuwname,'unwrap_ifg_index_sm');
      warning('on','MATLAB:load:variableNotFound');
      if isfield(phuw,'unwrap_ifg_index_sm');
        unwrap_ifg_index=phuw.unwrap_ifg_index_sm;
      else
        unwrap_ifg_index=[1:ps.n_image];
      end
    end

    if length(value_type)>2 & (value_type(1:3)=='usb'|value_type(1:3)=='rsb') & isempty(ifg_list)
        ifg_list=unwrap_ifg_index_sb;
    end
else
    if strcmpi(unwrap_ifg_index,'all') 
        unwrap_ifg_index=[1:ps.n_ifg];
    end
end
if (value_type(1)=='u' | value_type(1)=='a') & isempty(ifg_list)
    ifg_list=unwrap_ifg_index;
end

    


switch(value_type)
    case {'w'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc
     case {'w-d'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ~strcmpi(getparm('small_baseline_flag'),'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        ph_all=ph_all.*exp(-j*scla.ph_scla);
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
     case {'w-o'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ~strcmpi(getparm('small_baseline_flag'),'y')
            scla=load(sclaname);
        else
            scla=load(sclasbname);
        end
        ph_all=ph_all.*exp(-j*scla.ph_ramp);
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
     case {'w-do'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ~strcmpi(getparm('small_baseline_flag'),'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        ph_all=ph_all.*exp(-j*(scla.ph_scla+scla.ph_ramp));
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
     case {'w-dm'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        scla=load(sclasmoothname);
        ph_all=ph_all.*exp(-j*scla.ph_scla);
        ph_all=ph_all.*repmat(exp(-j*scla.C_ps_uw),1,size(ph_all,2));
        ph_all(:,ps.master_ix)=1;
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
     case {'w-dmo'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        scla=load(sclasmoothname);
        ph_all=ph_all.*exp(-j*(scla.ph_scla+scla.ph_ramp));
        ph_all=ph_all.*repmat(exp(-j*scla.C_ps_uw),1,size(ph_all,2));
        ph_all(:,ps.master_ix)=1;
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
    case {'u'}
        phuw=load(phuwname);
        ph_all=phuw.ph_uw;
        clear phuw
        ref_ps=ps_setref;
    case {'usb'}
        uw=load(phuwsbname);
        ph_all=uw.ph_uw;
        clear uw
        ref_ps=ps_setref;
        textsize=0;
    case {'dsb'}
        scla=load(sclasbname);
        ph_all=scla.K_ps_uw;
        clear scla
        ref_ps=ps_setref;
        units='rad/m';
    case {'usb-d'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        ph_all=uw.ph_uw - scla.ph_scla;
        clear uw scla
        ref_ps=ps_setref;
        textsize=0;
    case {'usb-o'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        ph_all=uw.ph_uw - scla.ph_ramp;
        clear uw scla
        ref_ps=ps_setref;
        textsize=0;
    case {'usb-do'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        ph_all=uw.ph_uw - scla.ph_scla - scla.ph_ramp;
        clear uw scla
        ref_ps=ps_setref;
        textsize=0;
    case {'usb-a'}
        uw=load(phuwsbname);
        aps=load(apssbname);
        ph_all=uw.ph_uw - aps.ph_aps_slave;
        clear uw aps
        ref_ps=ps_setref;
        textsize=0;
    case {'usb-da'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        aps=load(apssbname);
        ph_all=uw.ph_uw - scla.ph_scla - aps.ph_aps_slave;
        clear uw scla aps
        ref_ps=ps_setref;
        textsize=0;
    case {'rsb'}
        uw=load(phuwsbname);
        res=load(phuwsbresname);
        ph_all=zeros(size(uw.ph_uw));
        ph_all(:,unwrap_ifg_index_sb)=uw.ph_uw(:,unwrap_ifg_index_sb)-res.ph_res(:,unwrap_ifg_index_sb);
        clear uw
        ref_ps=ps_setref;
        textsize=0;
    case {'asb'}
        aps=load(apssbname);
        ph_all=aps.ph_aps_slave;
        clear aps
        ref_ps=ps_setref;
        textsize=0;
    case {'u-dms'}
        uw=load(phuwname);
        scn=load(scnname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scn.ph_scn_slave - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla;
        clear uw scn scla
        ph_all(:,ps.master_ix)=0;
        ref_ps=ps_setref;
    case {'u-dmas'}
        uw=load(phuwname);
        scn=load(scnname);
        scla=load(sclaname);
        aps=load(apsname);
        ph_all=uw.ph_uw - scn.ph_scn_slave - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla - aps.ph_aps_slave;
        clear uw scn scla
        ref_ps=ps_setref;
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
    case {'u-dmo'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla - scla.ph_ramp;
        ph_all(:,master_ix)=0;
        clear uw scla
        ref_ps=ps_setref;
    case {'u-dma'}
        uw=load(phuwname);
        scla=load(sclaname);
        aps=load(apsname);
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla - aps.ph_aps_slave;
        ph_all(:,master_ix)=0;
        clear uw scla aps
        ref_ps=ps_setref;
    case {'u-a'}
        uw=load(phuwname);
        aps=load(apsname);
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - aps.ph_aps_slave;
        ph_all(:,master_ix)=0;
        clear uw aps
        ref_ps=ps_setref;
    case {'u-d'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scla.ph_scla;
        clear uw scla
        ref_ps=ps_setref;
    case {'u-o'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scla.ph_ramp;
        clear uw scla
        ref_ps=ps_setref;
    case {'u-do'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scla.ph_scla - scla.ph_ramp;
        clear uw scla
        ref_ps=ps_setref;
    case {'u-m'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw;
        ph_all=uw.ph_uw - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2));
        ph_all(:,master_ix)=0;
        clear uw scla
        ref_ps=ps_setref;
    case {'a'}
        aps=load(apsname);
        %ph_all=exp(j*ph_scn);
        ph_all=aps.ph_aps_slave;
        clear aps
        ref_ps=ps_setref;
    case {'s'}
        scn=load(scnname);
        %ph_all=exp(j*ph_scn);
        ph_all=scn.ph_scn_slave;
        clear scn
        ref_ps=ps_setref;
    case {'m'}
        %scn=load(scnname);
        %ph_all=scn.ph_scn_master;
        %clear scn
        scla=load(sclaname);
        ph_all=scla.C_ps_uw;
        clear scla
        ref_ps=ps_setref;
    case {'d'}
        scla=load(sclaname);
        ph_all=scla.K_ps_uw;
        clear scla
        ref_ps=ps_setref;
        units='rad/m';
    case {'o'}
        scla=load(sclaname);
        ph_all=scla.ph_ramp;
        clear scla
        ref_ps=ps_setref;
        units='rad/m';
    case {'v'}
        uw=load(phuwname);
        ph_uw=uw.ph_uw;
        clear uw
        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        % Each ifg has master APS - slave APS, including master 
        % (where slave APS = master APS) so OK to include master in inversion
        G=[ones(size(day)),day-master_day]; 
        
        m=G\double(ph_uw');
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        textsize=0;
        units='mm/yr';
     case {'v-d'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_uw=uw.ph_uw - scla.ph_scla;
        clear uw scla
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        G=[ones(size(day)),day-master_day] ;
        
        m=G\double(ph_uw');
        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        %ph_all=m(3,:)'/4/pi*lambda*1000;
        textsize=0;
        units='mm/yr';
     case {'v-do'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_uw=uw.ph_uw - scla.ph_ramp - scla.ph_scla;
        clear uw scla
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        G=[ones(size(day)),day-master_day] ;
        
        m=G\double(ph_uw');
        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        textsize=0;
        units='mm/yr';
     case {'v-o'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_uw=uw.ph_uw - scla.ph_ramp;
        clear uw scla
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        G=[ones(size(day)),day-master_day] ;
        
        m=G\double(ph_uw');
        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        textsize=0;
        units='mm/yr';
     case {'v-da'}
        uw=load(phuwname);
        scla=load(sclaname);
        aps=load(apsname);
        ph_uw=uw.ph_uw - scla.ph_scla- aps.ph_aps_slave;
        clear uw scla aps
        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        G=[ones(size(day)),day-master_day]; 
        
        m=G\double(ph_uw');
        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        textsize=0;
        units='mm/yr';
     case {'v-a'}
        uw=load(phuwname);
        aps=load(apsname);
        ph_uw=uw.ph_uw - aps.ph_aps_slave;
        clear uw aps
        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        G=[ones(size(day)),day-master_day]; 
        
        m=G\double(ph_uw');
        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        %ph_all=m(3,:)'/4/pi*lambda*1000;
        textsize=0;
        units='mm/yr';
     case {'v-ds'}
        uw=load(phuwname);
        scla=load(sclaname);
        scn=load(scnname);
        ph_uw=uw.ph_uw - scla.ph_scla - scn.ph_scn_slave;
        clear uw scla scn

        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        G=[ones(size(day)),day-master_day]; 
        
        m=G\double(ph_uw');
        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        textsize=0;
        units='mm/yr';
     case {'v-das'}
        uw=load(phuwname);
        scla=load(sclaname);
        scn=load(scnname);
        aps=load(apsname);
        ph_uw=uw.ph_uw - scla.ph_scla- aps.ph_aps_slave - scn.ph_scn_slave;
        clear uw scla aps scn

        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        G=[ones(size(day)),day-master_day]; 
        
        m=G\double(ph_uw');
        try
            save mean_v m
        catch
            fprintf('Warning: Read access only, velocities were not saved\n')
        end
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        textsize=0;
        units='mm/yr';
   case {'vdrop-d'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_uw=uw.ph_uw - scla.ph_scla;
        clear uw scla
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index);
        
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        G=[ones(size(day)),day-master_day] ;
        n=size(ph_uw,2);
        lambda=getparm('lambda');
        ph_all=zeros(size(ph_uw));
        for i=1:n
            m=G([1:i-1,i+1:end],:)\double(ph_uw(:,[1:i-1,i+1:n])');
            ph_all(:,i)=-m(2,:)'*365.25/4/pi*lambda*1000; 
        end
        textsize=0;     
        units='mm/yr';
   case {'vsb'}
        phuw=load(phuwsbname);
        ref_ps=ps_setref;
        unwrap_ifg_index=getparm('unwrap_ifg_index')
        if strcmp(unwrap_ifg_index,'all')
            unwrap_ifg_index=[1:n_ifg];
        end	
        ph_uw=phuw.ph_uw(:,unwrap_ifg_index);
        ifgday_ix=ps.ifgday_ix(unwrap_ifg_index,:);
        clear phuw
        ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),n_ps,1);
        G=[ones(size(ifgday_ix(:,1))),day(ifgday_ix(:,2))-day(ifgday_ix(:,1))];
        m=G\double(ph_uw');
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; 
        textsize=0;
        units='mm/yr';
   case {'vsb-d'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        ph_all=uw.ph_uw - scla.ph_scla;
        clear uw scla
        ref_ps=ps_setref;
        unwrap_ifg_index=getparm('unwrap_ifg_index')
        if strcmp(unwrap_ifg_index,'all')
            unwrap_ifg_index=[1:n_ifg];
        end	
        ph_all=ph_all(:,unwrap_ifg_index);
        ifgday_ix=ps.ifgday_ix(unwrap_ifg_index,:);
        clear phuw
        ph_all=ph_all-repmat(mean(ph_all(ref_ps,:)),n_ps,1);
        G=[ones(size(ifgday_ix(:,1))),day(ifgday_ix(:,2))-day(ifgday_ix(:,1))];
        m=G\double(ph_all');
        lambda=getparm('lambda');
        ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; 
        textsize=0;
        units='mm/yr';
    case {'p'}
        pm=load(pmname);
        %pm=load(rcuwname);
        ph_all=pm.ph_patch./abs(pm.ph_patch);
        if n_ifg~=size(ph_all,2)
            ph_all=[ph_all(:,1:ps.master_ix-1),zeros(ps.n_ps,1),ph_all(:,ps.master_ix:end)];
        end
        clear pm
        if ref_ifg~=0
            ph_all=ph_all.*repmat(conj(ph_all(:,ref_ifg)),1,n_ifg);
            ph_all(:,ref_ifg)=1; % may not be so because of rounding errors
        end
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
    case {'vs-d'}
        %if exist(['./',meanvname,'.mat'],'file')==0
            ps_mean_v(ifg_list,200,0);
        %end
        ifg_list=[];
        mv=load(meanvname);
        ph_all=mv.mean_v_std;
        units='mm/yr';
    case {'vs-do'}
        ps_mean_v(ifg_list,200,1);
        ifg_list=[];
        mv=load(meanvname);
        ph_all=mv.mean_v_std;
        units='mm/yr';

        
        
    otherwise
        error('unknown value type')
end

if isempty(ifg_list)
    ifg_list=1:size(ph_all,2);
end
n_ifg_plot=length(ifg_list);


[Y,X]=meshgrid([0.7:-0.2:0.1],[0.1:0.1:0.8]);
xgap=0.1;
ygap=0.2;

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
else
    n_y=ceil((n_ifg_plot)/n_x);
end


d_x=useratio/n_x;
d_y=d_x/ar*fig_ar;
if d_y>useratio/n_y
    d_y=useratio/n_y; 
    d_x=d_y*ar/fig_ar;
end

h_y=0.95*d_y;
h_x=h_y*ar/fig_ar;
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
            mean_ph(i)=mean(ref_ph(~isnan(ref_ph(:,i)),i));
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
  figure
  set(gcf,'renderer','zbuffer')
  i_im=0;
  for i=ifg_list
    %subplot(5,7,i);
    i_im=i_im+1;
    if n_ifg_plot>1
        axes('position',[imX(i_im),imY(i_im),h_x,h_y])
    end

    %axes('position',[X(i),Y(i),xgap,ygap])
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
        else
            h=colorbar('SouthOutside');
        end
        xlim=get(h,'xlim');
        set(h,'xlim',[xlim(2)-64,xlim(2)])
        if diff(lims)>1 | diff(lims)==0
            plotlims=round(lims*10)/10;
        else
            limorder=ceil(-log10(diff(lims)))+2;
            plotlims=round(lims*10^limorder)/10^limorder;
        end
        set(h,'xtick',[xlim(2)-64,xlim(2)],'Xticklabel',plotlims,'xcolor',textcolor,'ycolor',textcolor,'fontweight','bold','color',textcolor)
        h=xlabel(h,units);
        pos=get(h,'position');
        pos(2)=pos(2)/2.2;
        set(h,'position',pos);
    end

  end
end
    

fprintf('Color Range: %g to %g %s\n',lims,units)


    
    

