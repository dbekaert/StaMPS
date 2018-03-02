function [h_fig,lims,ifg_data_RMSE,h_axes_all]=ps_plot(value_type,varargin)
%function [h_fig,lims,ifg_data_RMSE,h_axes_all]=ps_plot(value_type,plot_flag,lims,ref_ifg,ifg_list,n_x,...
%    cbar_flag,textsize,textcolor,lon_rg,lat_rg)
% PS_PLOT plot ps values for selected ifgs
%    PS_PLOT(VALUE_TYPE,BACKGROUND,PHASE_LIMS,REF_IFG,IFG_LIST,N_X,...
%            CBAR_FLAG,TEXTSIZE,TEXTCOLOR,LON_RG,LAT_RG,UNITS) 
%
%    In the case of phase, +ve values imply displacement away from the satellite
%       when the master is earlier than the slave.
%    In the case of velocities, +ve values are towards the satellite.
%
%    an invalid VALUE_TYPE related to double removing of atmosphere:
%    anything when removing master and estimated tropospheric signals: 'ma'
%
%    valid VALUE_TYPE's are:
%%%% BASIC
%    'hgt'  for topography
%    'w'    for wrapped phase
%    'p'    for spatially filtered wrapped phase 
%    'u'    for unwrapped phase
%    'usb'  for unwrapped phase of small baseline ifgs
%    'rsb'  for residual between unwrapped phase of sb ifgs and inverted
%    'd'    for spatially correlated DEM error (rad/m) estimated from SM or from SB depending on your processing
%    'D'    for spatially correlated DEM error (rad/m) from SM inverted data when using SB
%    'm'    for AOE phase due to master
%    'o'    for orbital ramps 
%    's'    for atmosphere and orbit error (AOE) phase due to slave (not compatible with 'a')
%    'a'*   for stratisfied topo-correlated atmosphere using TRAIN dependency. 
%    'asb'* for stratisfied topo-correlated atmosphere of small baselines using TRAIN dependency
%    *'a' & 'asb' require an additional string to indicate which APS correction to be applied:
%             -'a_linear'                                 = (mixed) Tropospheric aps correction using linear phase-topo correction
%             -'a_powelaw'                                = (mixed) Tropospheric aps correction using powerlaw phase-topo relationship
%             -'a_meris'                                  = (wet) tropospheric aps correction using MERIS data (wet)
%             -'a_modis'                                  = (wet) tropospheric aps correction using MODIS data (wet)
%             -'a_modisRecal'                             = (wet) tropospheric aps correction using MODIS recalibrated data
%             -'a_erai','a_erai-wet','a_erai-hydro'       = (total, wet, hydro) tropospheric aps correction using ERA-I data
%             -'a_wrf','a_wrf-wet','a_wrf-hydro'          = (total, wet, hydro) tropospheric aps correction using WRF model data
%             -'a_merra','a_merra-wet','a_merra-hydro'    = (total, wet, hydro) tropospheric aps correction using MERRA model data
%             -'a_merra2','a_merra2-wet','a_merra2-hydro' = (total, wet, hydro) tropospheric aps correction using MERRA-2 model data
%             -'a_gacos'                                  = (total) tropospheric aps correction using GACOS data
%             -'a_pk'                                     = spatial transfer coefficient K for interferograms or specific bands.
%    'v'    for mean LOS velocity (MLV) in mm/yr estimated from SM or from SB depending on your processing
%    'V'    for mean LOS velocity  of MLV (mm/yr) estimated from SM inverted data when using SB
%    'vs'   for standard deviation of MLV (mm/yr) estimated from SM or from SB depending on your processing
%    'Vs'   for standard deviation of MLV (mm/yr) estimated from SM inverted data when using SB
%    'vdrop'for  MLV calculated from all but current ifg (mm/yr)

%%%% CORRECTIONS to BASIC
%    'w-a'  for wrapped phase corrected for atmopshere
%    'w-d'  for wrapped phase minus smoothed dem error
%    'w-o'  for wrapped phase minus orbital ramps
%           also 'w-dm','w-do','w-dmo'
%    'u-d'  for unwrapped phase minus dem error  
%    'u-m'  for unwrapped phase minus and master AOE  
%    'u-o'  for unwrapped phase minus orbital ramps 
%    'u-a'  for unwrapped phase minus stratisfied topo-correlated atmosphere
%           also 'u-dm','u-do','u-da','u-ao','u-dmo','u-dma','u-dms','u-dmao','u-dmos'
%    'usb-d' 
%    'usb-o' 
%    'usb-a' 
%           also 'usb-do' ,'usb-da','usb-dao'
%    'v-d'  for mean LOS velocity (MLV) in mm/yr minus smoothed dem error
%    'v-o' 
%    'v-a' 
%           also 'v-do','v-da',v-dao','v-s','v-so'and the capital V options
%    'vs-d' for standard deviation of MLV (mm/yr) minus smoothed dem error
%    'vs-o' 
%    'vs-a' 
%      also 'vs-do', 'vs-da', 'vs-dao' and the capital Vs options
%    'vdrop-d' 
%    'vdrop-o' 
%      also 'vdrop-do' 
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
%             Note needs to be SM (v-option) of SB inverted to SM
%             time-series (V-option)
%
%
%    'ifg i' = only for 'a_p' show the topopgraphy correlated aps correction
%              for the ith interferogram for all spatial bands. Note that
%              the definition of ifg_list changes the spatial bands for this option.
%
%   ADDITONAL DATA DISPLAY
%     Plot in additional also other geocoded information, e.g stations or 
%     LOS data on top of the interferograms. To use this functionality use 
%     the flag:  'ext PATH', with PATH the full path to the data location.
%     To display identical information in each interferogram include a mat file
%     on this path. In case the to be displayed data varies depending on
%     the ifg dates, then only specify the path. At this path the data should 
%     be storred in date1_date2.mat files, where the dates are in YYYYMMDD format. 
%     The file itself should contain at least a lonlat variable. When also a 
%     ph_disp (phase) variable is given the markers will be collored acordingly. 
%     NOTE: - Not each interferogram date needs to have a file asociated with it.
%           - External ph_disp data are fixed to the ifgs by minimizing the mean residual.
%           - In case a reference area is selected, it needs to have data coverage in all datasets.
%           - Markers are given as by default by collored squares. In case
%             they were saturated by the colorbar they become circular
%             markers.
%
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
%   12/2010 KS: Added check to see if orbital ramp is calculated in -o
%   12/2010 AH: fix vdrop
%   03/2011 AH: Change ts plot radius to m
%   03/2011 AH: Save files to home directory if read only directory
%   03/2011 DB: Remove raster lines when plotting time series
%   10/2011 AH: Check that pm file exists before attempting to load
%   01/2012 AH: Subtract master AOE for 'ts' plot
%   01/2012 AH: Remove code to subtract SULA error from 'd' plots
%   12/2012 AH: plot raw phase if psver=1
%   01/2013 DB: Fix plotting of data in case of SB directory
%   01/2013 AH: Add topo-correlated atmosphere options
%   04/2013 DB: Add figure handle as output
%   04/2013 DB: Add colorbar axis extremes as output
%   04/2013 DB: Include plotting of atmospheric correction from meris data
%   04/2013 DB: Topo correlated aps correction using 'a_m','a_l','a_p' options 
%   04/2013 DB: Adding the same velocity plot options for SB as for SM
%   04/2013 DB: Bug fix for aps, include option to show aps_p for all bands
%               for infividual interferograms, by including 'ifg i' option
%   04/2013 DB: Include an option to show in additional also other data
%               like e.g. LOS GPS displacements
%   04/2013 DB: Modify such that for bandfiltered data the external data
%               for that date alone is shown.
%   04/2013 DB: Fix the plotting of external data by minimizing the
%               residual with the insar. Output the RMSE for each ifgs.
%   05/2013 DB: Plot the ifg number for SB ifgs
%   05/2013 DB: Allow inputdata to be saved.
%   05/2013 DB: Option to plot topography
%   05/2013 DB: When full filepath is given, plot same ext data for each ifgs
%               and output the RMSE information.
%   06/2013 DB: Allow units to be specified and fix wrong unit for 'o' option
%   09/2013 DB: Added extra options and fix warning for bandfilter option
%   10/2013 DB: Included ERA-I tca option
%   11/2013 DB: Incorporate the inversion of the troposphere from SB to PS
%   11/2013 DB: Including WRF atmospheric option
%   03/2014 AH: vs-a etc added
%   04/2014 DB: Subtract envisat oscialtor drift on the go for all emvisat
%               plotting options
%   05/2014 DB: When saving give empty outputs 
%   05/2014 DB: Fix oscialtor drift when plotting 'u' and 'usb' options
%   05/2014 DB: Add extra output argument h_axes_all
%   05/2014 DB: Hardcode topogrpahy to have a GMT colormap for relief
%   06/2014 DB: Fix bug in case of plotting u
%   07/2014 DB: Fix a bug, which already corrected itself latter on
%   07/2014 DB: Fix in case there is no data in the reference area
%   08/2014 DB: Bug fix for ps_plot V-DOA
%   08/2014 DB: Allow to plot w-a, and wa (wrapped atmosphere)
%   08/2014 DB: Include modis recalibrated data support
%   10/2014 DB: Include support for ionopsheric delays
%   03/2015 DB: Remove the reference for the K spatial map option
%   03/2015 DB: Update and clean for relaese with TRAIN
%   05/2015 DB: Fix to separate hgt from t options
%   12/2016 DB: Update output syntax, remove warning ps_plot('V-Do')
%   06/2016 DB: Fix in case the master is still included in SM inversion
%   11/2017 DB: Updating the TRAIN plotting conventions and adding new options
%               Only let ts option run through for option ('v') and ('V')
%   11/2017 DB: make the s option work again
%   ======================================================================

stdargin = nargin ; 
parseplotprm  % check if 'ts', 'a_m', 'a_l', 'a_e' ('a_eh' for hydrostatic and 'a_ew' for wet), 'a_w' ('a_wh' for hydrostatic and 'a_ww' for wet),
% 'a_p' ('a_pk' for spatial map of coefficent K), 'i_as', 'ifg i^th', 'ext PATH ' is specified

reference_flag=0;

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
    if plot_flag==1 | plot_flag==2 |plot_flag==6
        textcolor=[0 0 0.004];
        textcolor2=textcolor;

    else
        textcolor=[1 1 0.996];
        textcolor2=[0 0 0.004];
    end
else
    textcolor2=textcolor;
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

if stdargin < 12
    units=[];
end
% reference radius in case of external data to be plotted
ref_radius_data=1000;

n_y=0;

load psver
psname=['./ps',num2str(psver)];
pmname=['./pm',num2str(psver)];
rcname=['./rc',num2str(psver),'.mat'];
phname=['./ph',num2str(psver),'.mat'];
phuwname=['./phuw',num2str(psver)];
phuwsbname=['./phuw_sb',num2str(psver)];
phuwsbresname=['./phuw_sb_res',num2str(psver)];
scnname=['./scn',num2str(psver)];
ifgstdname=['./ifgstd',num2str(psver)];
apsbandsname = ['./tca_bands' num2str(psver) '.mat'];
apsbandssbname = ['./tca_bands_sb' num2str(psver) '.mat'];
apsname=['./tca',num2str(psver)];
apssbname=['./tca_sb',num2str(psver)];
iononame = ['./ica' num2str(psver) '.mat'];
ionosbname = ['./ica_sb' num2str(psver) '.mat'];
tidename=['./tide',num2str(psver)];
tidesbname=['./tide_sb',num2str(psver)];
hgtname=['./hgt',num2str(psver)];


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
scla_deramp = getparm('scla_deramp');


% Making dummy files with the data when displaying band filtered data
if aps_band_flag==1 
    % do not use the parm list to reject images as it is meant for ifgs.
    drop_ifg_index = [];
    % put the colorbar at the first chosen frequency band
    if isempty(ifg_list)~=1
        master_ix = ifg_list(1);
    else
        master_ix = 1;
    end
    % getting the spatial extends of the bands
    if exist('parms_aps.mat')==2
        bands = getparm_aps('powerlaw_spatial_bands');
    else
        bands = [];
    end
    
   
    % the data        
    if strcmp(small_baseline_flag,'y')         
        % keep only the ith ifg for aps
        aps_bands = load(apsbandssbname);
        apssbname=['./tca_sb_temp',num2str(psver)];

        if aps_flag ==11
            eval(['K_tropo_powerlaw =  aps_bands.K_tropo_powerlaw_ifg_' num2str(ifg_number) ';']);            
            save(apssbname,'K_tropo_powerlaw');
            n_bands = size(K_tropo_powerlaw,2); 
        else
            eval(['ph_tropo_powerlaw =  aps_bands.ph_tropo_powerlaw_ifg_' num2str(ifg_number) ';']);
            save(apssbname,'ph_tropo_powerlaw');
            n_bands = size(ph_tropo_powerlaw,2); 
        end        
        delete_temp_files{1} = apssbname;
            
        % for phase        
        counter = 2;
        if strcmp(value_type(1:3),'usb')
            % loading the phase
            uw=load(phuwsbname);
            ph_all=uw.ph_uw;
            phuwsbname=['./phuw_sb_temp',num2str(psver)];
            % selecting the interferogram for display
            ph_uw = repmat(ph_all(:,ifg_number),1,n_bands);
            save(phuwsbname,'ph_uw');
            clear ph_uw   
            delete_temp_files{counter} = phuwsbname;
            counter = counter+1;
        end

        
        if strcmp(value_type,'usb') || strcmp(value_type,'asb')
        else
            % other options still need to be included
            error('This option is not included in ps_plot.m')
        end
        
       
        
    else
        % keep only the ith ifg for aps
        aps_bands = load(apsbandsname);
        apsname=['./tca_temp',num2str(psver)];
        
        if aps_flag ==11
            eval(['K_tropo_powerlaw =  aps_bands.K_tropo_powerlaw_ifg_' num2str(ifg_number) ';']);            
            save(apsname,'K_tropo_powerlaw');           
            n_bands = size(K_tropo_powerlaw,2); 
        else
            eval(['ph_tropo_powerlaw =  aps_bands.ph_tropo_powerlaw_ifg_' num2str(ifg_number) ';']);
            save(apsname,'ph_tropo_powerlaw');
            n_bands = size(ph_tropo_powerlaw,2); 
        end
        delete_temp_files{1} = apsname;

        % for phase        
        counter = 2;
        if strcmp(value_type(1),'u')
            % loading the phase
            uw=load(phuwname);
            ph_all=uw.ph_uw;
            phuwname=['./phuw_temp',num2str(psver)];
            % selecting the interferogram for display
            ph_uw = repmat(ph_all(:,ifg_number),1,n_bands);
            save(phuwname,'ph_uw');
            clear ph_uw   
            delete_temp_files{counter} = phuwname;
            counter = counter+1;
        end
        
        if strcmp(value_type,'u') || strcmp(value_type,'a')
        else
            % other options still need to be included
            error('This option is not included in ps_plot.m')
        end
        
        
        
        
    end
   
    fig_name_suffix = [' bandfiltered data of ifg ' num2str(ifg_number)];

else
    fig_name_suffix = '';
    bands = [];
end



if aps_band_flag==1 & isempty(ifg_list)
    % Change the definition of ifg_list from ifgs to the different bands
    ifg_list = [1:n_bands];
else
    if strcmpi(small_baseline_flag,'y')
        unwrap_ifg_index_sb=setdiff([1:ps.n_ifg],drop_ifg_index);
        if ischar(value_type) & value_type(1)~='w' & value_type(1)~='p' & exist([phuwname,'.mat'],'file')
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
        if ischar(value_type)~=1 & isempty(ifg_list) % [DB] to plot data matrix which is not a valuetype
                    ifg_list=unwrap_ifg_index_sb;
        elseif ischar(value_type)==1 & length(value_type)>2 & (value_type(1:3)=='usb'|value_type(1:3)=='rsb'| value_type(1:3)=='asb' | value_type(1:3)=='isb') & isempty(ifg_list)
                    ifg_list=unwrap_ifg_index_sb;
        end
    else
        unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);
    end
    
    if (value_type(1)=='u' | value_type(1)=='s' | value_type(1)=='a' | value_type(1)=='i' | value_type(1)=='w') & isempty(ifg_list)
        ifg_list=unwrap_ifg_index;
    end
    if ischar(value_type)~=1 & size(value_type,2)<length(ifg_list)
      ifg_list = [1:size(value_type,2)]';
    end 
end

if isempty(units) && ischar(value_type)==1
      units='rad'; 
end

% flag for oscilator drift
forced_sm_flag=0;
if ischar(value_type)==1
group_type=value_type;
if length(group_type)>1 & strcmpi(group_type(1:2),'vs')
    group_type='vs'; 
    if strcmpi(small_baseline_flag,'y') & strcmp(value_type(1:2),'vs')
        use_small_baselines=1;
        fprintf('Velocity std devs calculated from small baseline interferograms\n')
    else
        use_small_baselines=0;
        forced_sm_flag=1;
    end
elseif strcmpi(group_type(1),'v')
    if strcmpi(small_baseline_flag,'y') & strcmp(group_type(1),'v')
        group_type='vsb';
        fprintf('Velocities calculated from small baseline interferograms\n')
    else
        group_type='v';
        %TEMPsb_invert_uw;               % Added to make sure if the network is changed it is updated
        forced_sm_flag=1;
    end
elseif strcmpi(group_type(1),'d')
    if strcmpi(small_baseline_flag,'y') & strcmp(group_type(1),'d')
        group_type='dsb';
        fprintf('DEM error calculated from small baseline interferograms\n')
    else
        group_type='d';
        forced_sm_flag=1;
    end
end

% only allow TS option for 'v' and 'V' options
if ts_flag==1 && ~(strcmp(group_type,'v'))
    fprintf('\n\n''ts'' option can only be called for:\n    -SM: ''v''-options\n    -SB to SM inverted: ''V''-options\n')
    error('Incorrect value-type for ''ts'' option')
end


if strcmp(value_type(1),'u') || strcmp(value_type(1),'a')
   if length(value_type)>=3
       if strcmp(value_type(1:3),'usb') || strcmp(value_type(1:3),'asb')
           forced_sm_flag=0;
       else
           forced_sm_flag=1;
       end
   else
       forced_sm_flag=1;       
   end

end


% related to the oscilator drift correction - envisat only
if forced_sm_flag==1
    [ph_unw_eni_osci,v_envi_osci] = env_oscilator_corr([],forced_sm_flag);
else
    % Compute the envisat oscialator drift when needed
    [ph_unw_eni_osci,v_envi_osci] = env_oscilator_corr;
end

% Separate tide t from topography hgt option
value_type=lower(value_type);
aps_tide_iono_flag = 1;             % do the checks anyway, but make sure its not hgt option
if length(value_type)>=3
    if strcmpi(value_type(1:3),'hgt')
        aps_tide_iono_flag = 0;
    end
end
if aps_tide_iono_flag==1
    % Check if the tropopsheric needs to be inverted from SB to PS
    if strcmpi(small_baseline_flag,'y') && ~isempty(strfind(value_type,'a')) && isempty(strfind(value_type,'sb')) &&  (strcmpi(group_type,'vsb'))==0
        if isempty(strfind(value_type,'w')) 
            if (strcmpi(group_type,'vs'))==1 && use_small_baselines==0
               % this is small baselines velocity with SB aps correction
               %TEMP
               sb_invert_aps(aps_flag);
            elseif (strcmpi(group_type,'vs'))~=1
               % this is small baselines velocity with SB aps correction
               %TEMP
               sb_invert_aps(aps_flag);
            end
        end
    end

    % Check if the tidal correction needs to be inverted from SB to PS
    if strcmpi(small_baseline_flag,'y') && ~isempty(strfind(value_type,'t')) && isempty(strfind(value_type,'sb')) &&  (strcmpi(group_type,'vsb'))==0     
        if isempty(strfind(value_type,'w')) 
            if (strcmpi(group_type,'vs'))==1 && use_small_baselines==0
               % this is small baselines velocity with SB tide correction
               sb_invert_tide;
            elseif (strcmpi(group_type,'vs'))~=1
               % this is small baselines velocity with SB tide correction
               sb_invert_tide;
            end
        end
    end
    
    % Check if the ionopshere need to be inverted from SB to PS
    if strcmpi(small_baseline_flag,'y') && ~isempty(strfind(value_type,'i')) && isempty(strfind(value_type,'sb')) &&  (strcmpi(group_type,'vsb'))==0
        if isempty(strfind(value_type,'w')) 
            if (strcmpi(group_type,'vs'))==1 && use_small_baselines==0
               % this is small baselines velocity with SB aps correction
               sb_invert_iono(iono_flag);
            elseif (strcmpi(group_type,'vs'))~=1
               % this is small baselines velocity with SB aps correction
               sb_invert_iono(iono_flag);
            end
        end
    end
end

switch(group_type)
    case {'hgt'}
        hgt = load(hgtname);
        ph_all=hgt.hgt;
        clear rc
        fig_name = 'hgt';  
        units='m';
        
        plot_color_scheme_old = getparm('plot_color_scheme');
        setparm('plot_color_scheme','GMT_globe');
    case {'w'}
        if exist(rcname,'file')
            rc=load(rcname);
            ph_all=rc.ph_rc;
        else
            rc=load(phname);
            ph_all=rc.ph;
        end
        
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all.*exp(-j*ph_unw_eni_osci);
        
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc
        fig_name = 'w';
    case {'w-a'}
        if exist(rcname,'file')
            rc=load(rcname);
            ph_all=rc.ph_rc;
        else
            rc=load(phname);
            ph_all=rc.ph;
        end
        
        
        if strcmpi(small_baseline_flag,'y')
            aps=load(apssbname);
        else
            aps=load(apsname);
        end
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
       
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all.*exp(-j*ph_unw_eni_osci).*exp(-j*aps_corr);
        
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc
        clear  aps aps_corr
        fig_name = ['w-a' fig_name_tca];

     case {'w-da'}
        if exist(rcname,'file')
            rc=load(rcname);
            ph_all=rc.ph_rc;
        else
            rc=load(phname);
            ph_all=rc.ph;
        end
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        
        if strcmpi(small_baseline_flag,'y')
            aps=load(apssbname);
        else
            aps=load(apsname);
        end
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
       
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all.*exp(-j*ph_unw_eni_osci).*exp(-j*aps_corr).*exp(-j*scla.ph_scla);
        
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc
        clear  aps aps_corr
        fig_name = ['w-da' fig_name_tca];

        
         
        
    case {'wa'}       
        if strcmpi(small_baseline_flag,'y')
            aps=load(apssbname);
        else
            aps=load(apsname);
        end
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
       
        % subtract of oscialtor drift in case of envisat
        ph_all = exp(-j*aps_corr);
        ix = sum(ph_all~=1)==0;
       
        ph_all(:,ix)=complex(1,1);
        
        
        
        
        if ref_ifg~=0
%             ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear  aps aps_corr
        fig_name = ['wa' fig_name_tca];

        
     case {'w-d'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        ph_all=ph_all.*exp(-j*scla.ph_scla);
        
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all.*exp(-j*ph_unw_eni_osci);
        
        
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
        fig_name = 'w-d';
     case {'w-o'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        
        if strcmp('n',scla_deramp)
            disp('Warning: scla_deramp flag set to n. Set to y and rerun Step 7 before using the -o plot command.')
            return;
        end
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclaname);
        else
            scla=load(sclasbname);
        end
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all.*exp(-j*ph_unw_eni_osci);
        if sum(sum(ph_unw_eni_osci))~=0
            fprintf('Warning: note that the oscilator drift is also removed, make sure that the ramp is estimated after correction \n')
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
        if strcmp('n',scla_deramp)
            disp('Warning: scla_deramp flag set to n. Set to y and rerun Step 7 before using the -o plot command.')
            return;
        end
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        ph_all=ph_all.*exp(-j*(scla.ph_scla+scla.ph_ramp));
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all.*exp(-j*ph_unw_eni_osci);
        if sum(sum(ph_unw_eni_osci))~=0
            fprintf('Warning: note that the oscilator drift is also removed, make sure that the ramp is estimated after correction \n')
        end
        
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
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all.*exp(-j*ph_unw_eni_osci);
        
        ph_all(:,ps.master_ix)=1;
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
        fig_name = 'w-dm';
     case {'w-dmo'}
        rc=load(rcname);
        ph_all=rc.ph_rc;
        if strcmp('n',scla_deramp)
            disp('Warning: scla_deramp flag set to n. Set to y and rerun Step 7 before using the -o plot command.')
            return;
        end
        if ~strcmpi(small_baseline_flag,'y')
            scla=load(sclasmoothname);
        else
            scla=load(sclasbsmoothname);
        end
        ph_all=ph_all.*exp(-j*(scla.ph_scla+scla.ph_ramp));
        ph_all=ph_all.*repmat(exp(-j*scla.C_ps_uw),1,size(ph_all,2));
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all.*exp(-j*ph_unw_eni_osci);
        if sum(sum(ph_unw_eni_osci))~=0
            fprintf('Warning: note that the oscilator drift is also removed, make sure that the ramp is estimated after correction \n')
        end
        
        ph_all(:,ps.master_ix)=1;
        if ref_ifg~=0
            ph_all=ph_all.*conj(repmat(rc.ph_reref(:,ref_ifg),1,n_ifg));
        end
        clear rc scla
        fig_name = 'w-dmo';
    case {'u'}
        phuw=load(phuwname);
        ph_all=phuw.ph_uw;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;        
        clear phuw
        ref_ps=ps_setref;
        fig_name = 'u';
        units='rad';
         
    case {'usb'}
        uw=load(phuwsbname);
        ph_all=uw.ph_uw;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        clear uw
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb';
                        
        
    case {'dsb'}
        scla=load(sclasbname,'K_ps_uw');
        ph_all=scla.K_ps_uw;
        clear scla
        ref_ps=ps_setref;
        units='rad/m';
        fig_name = 'dsb';
    case {'usb-d'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        ph_all=uw.ph_uw - scla.ph_scla;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        clear uw scla
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb-d';
    case {'usb-o'}
        uw=load(phuwsbname);
        ph_all=uw.ph_uw;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
        clear uw scla
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb-o';
        units='rad';
    case {'usb-do'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        ph_all=uw.ph_uw - scla.ph_scla;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
        clear uw scla
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'usb-do';
    case {'usb-a'}
        uw=load(phuwsbname);
        aps=load(apssbname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['usb-a' fig_name_tca];
        ph_all=uw.ph_uw - aps_corr;
        clear uw aps aps_corr
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ref_ps=ps_setref;
        textsize=0;
    case {'usb-ai'}
        uw=load(phuwsbname);
        aps=load(apssbname);
        % removing troposphere
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        iono=load(ionosbname);
        % removing ionopshere
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        fig_name = ['usb-ai' fig_name_tca ' ' fig_name_ica];
        ph_all=uw.ph_uw - aps_corr - iono_corr;
        clear uw aps aps_corr
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        
        ref_ps=ps_setref;
        textsize=0;
     case {'usb-ao'}
        uw=load(phuwsbname);
        aps=load(apssbname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['usb-ao' fig_name_tca];
        ph_all=uw.ph_uw - aps_corr;
        clear uw aps aps_corr
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
        ref_ps=ps_setref;
        textsize=0;
    case {'usb-da'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        aps=load(apssbname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['usb-da' fig_name_tca];
        ph_all=uw.ph_uw - scla.ph_scla - aps_corr;
        clear uw scla aps aps_corr
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ref_ps=ps_setref;
        textsize=0;
        
     case {'usb-dai'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        aps=load(apssbname);
        iono=load(ionosbname);
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['usb-dai' fig_name_tca ' ' fig_name_ica] ;
        ph_all=uw.ph_uw - scla.ph_scla - aps_corr -iono_corr;
        clear uw scla aps aps_corr
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ref_ps=ps_setref;
        textsize=0;    
    case {'usb-dao'}
        uw=load(phuwsbname);
        scla=load(sclasbname);
        aps=load(apssbname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['usb-dao' fig_name_tca];
        ph_all=uw.ph_uw - scla.ph_scla - aps_corr;
        clear uw scla aps aps_corr
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
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
        fig_name = 'rsb';  
    case {'rsb-o'}
        uw=load(phuwsbname);
        res=load(phuwsbresname);
        ph_all=zeros(size(uw.ph_uw));
        ph_all(:,unwrap_ifg_index_sb)=uw.ph_uw(:,unwrap_ifg_index_sb)-res.ph_res(:,unwrap_ifg_index_sb);
        [ph_all] = ps_deramp(ps,ph_all );
        clear uw
        ref_ps=ps_setref;
        textsize=0;
        fig_name = 'rsb-o'; 
        
    case {'asb'}
        aps=load(apssbname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);    
        fig_name = ['asb' fig_name_tca];
        fig_name = [fig_name fig_name_suffix];
        ph_all=aps_corr;
        % when plotting spatial maps of K do not re-reference the data!
        if aps_flag~=11
            ref_ps=ps_setref;
        end
        textsize=0;
        
     case {'tsb'}
            tide=load(tidesbname);
            fig_name = ['tsb'];
            ph_all= tide.ph_tide;  
            ref_ps=ps_setref;
    case {'isb'}
        iono=load(ionosbname);
        %ph_all=exp(j*ph_scn);
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        fig_name = ['isb' fig_name_ica];
        fig_name = [fig_name fig_name_suffix];
        ph_all=iono_corr;
        clear iono iono_corr
        ref_ps=ps_setref;   
     case {'asb-o'}
        aps=load(apssbname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['asb-o' fig_name_tca];
        fig_name = [fig_name fig_name_suffix];
        ph_all=aps_corr;
        clear aps aps_corr
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
        ref_ps=ps_setref;
        textsize=0;       
    case {'u-dms'}
        uw=load(phuwname);
        scn=load(scnname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scn.ph_scn_slave - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla;
        clear uw scn scla
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ph_all(:,ps.master_ix)=0;
        ref_ps=ps_setref;
        fig_name = 'u-dms';
    case {'u-dm'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw; 
        %ph_all(:,unwrap_ifg_index)=ph_all(:,unwrap_ifg_index) - repmat(scla.C_ps_uw,1,length(unwrap_ifg_index));
        %ph_all(:,unwrap_ifg_index) = ph_all(:,unwrap_ifg_index) - scla.ph_scla(:,unwrap_ifg_index);
        ph_all=uw.ph_uw - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ph_all(:,master_ix)=0;
        clear uw scla
        ref_ps=ps_setref;
        fig_name = 'u-dm';
    case {'u-dmo'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla;
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
        ph_all(:,master_ix)=0;
        clear uw scla
        ref_ps=ps_setref;
        fig_name = 'u-dmo';
    case {'u-da'}
        uw=load(phuwname);
        scla=load(sclaname);
        aps=load(apsname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['u-da' fig_name_tca];
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - scla.ph_scla - aps_corr;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ph_all(:,master_ix)=0;
        clear uw scla aps aps_corr
        ref_ps=ps_setref;
        
        
      case {'u-dai'}
        uw=load(phuwname);
        scla=load(sclaname);
        aps=load(apsname);
        iono=load(iononame);
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['u-dai' fig_name_tca ' ' fig_name_ica];
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - scla.ph_scla - aps_corr - iono_corr;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ph_all(:,master_ix)=0;
        clear uw scla aps aps_corr
        ref_ps=ps_setref;
        
        
      case {'u-dait'}
        uw=load(phuwname);
        scla=load(sclaname);
        aps=load(apsname);
        iono=load(iononame);
        tide=load(tidename);
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['u-dait' fig_name_tca ' ' fig_name_ica];
        ph_all=uw.ph_uw; 
        ph_all=uw.ph_uw - scla.ph_scla - aps_corr - iono_corr- tide.ph_tide;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ph_all(:,master_ix)=0;
        clear uw scla aps aps_corr
        ref_ps=ps_setref;
 
    case {'u-dma'}
        error('This is not a valid option: you cannot remove atmosphere twice \n ')
    case {'u-dmao'}
        error('This is not a valid option: you cannot remove atmosphere twice \n ')

    case {'u-dmos'}
        uw=load(phuwname);
        scn=load(scnname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scn.ph_scn_slave - repmat(scla.C_ps_uw,1,size(uw.ph_uw,2)) - scla.ph_scla;
        clear uw scn scla
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
        ph_all(:,ps.master_ix)=0;
        ref_ps=ps_setref;
        fig_name = 'u-dms';
    case {'u-a'}
        uw=load(phuwname);
        aps=load(apsname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['u-a' fig_name_tca];
        ph_all=uw.ph_uw;
        ph_all=uw.ph_uw - aps_corr;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        if aps_band_flag~=1
            ph_all(:,master_ix)=0;
        end
        clear uw aps aps_corr
        ref_ps=ps_setref;
     case {'u-ai'}
        uw=load(phuwname);
        aps=load(apsname);
        iono=load(iononame);
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['u-ai' fig_name_tca ' ' fig_name_ica];
        ph_all=uw.ph_uw;
        ph_all=uw.ph_uw - aps_corr - iono_corr;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        if aps_band_flag~=1
            ph_all(:,master_ix)=0;
        end
        clear uw aps aps_corr
        ref_ps=ps_setref;  
        
        
     case {'u-ao'}
        uw=load(phuwname);
        aps=load(apsname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['u-ao' fig_name_tca];
        ph_all=uw.ph_uw;
        ph_all=uw.ph_uw - aps_corr;
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        if aps_band_flag~=1
            ph_all(:,master_ix)=0;
        end
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
        clear uw aps aps_corr
        ref_ps=ps_setref;
        
    case {'u-dao'}
        uw=load(phuwname);
        aps=load(apsname);
        scla=load(sclaname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['u-dao' fig_name_tca];
        ph_all=uw.ph_uw;
        ph_all=uw.ph_uw - aps_corr - scla.ph_scla;
        if aps_band_flag~=1
            ph_all(:,master_ix)=0;
        end
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
        clear uw aps aps_corr scla
        ref_ps=ps_setref;     
        
    case {'u-d'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all=uw.ph_uw - scla.ph_scla;
        clear uw scla
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ref_ps=ps_setref;
        fig_name = 'u-d';
        
    case {'u-di'}
        uw=load(phuwname);
        scla=load(sclaname);
        iono=load(iononame);
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        fig_name = ['u-di' fig_name_ica];
        fig_name = [fig_name fig_name_suffix];
        ph_all=uw.ph_uw - scla.ph_scla - iono_corr;
        clear uw scla
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ref_ps=ps_setref;
        
    case {'u-dit'}
        uw=load(phuwname);
        scla=load(sclaname);
        iono=load(iononame);
        tide=load(tidename);
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        fig_name = ['u-dit' fig_name_ica];
        fig_name = [fig_name fig_name_suffix];
        ph_all=uw.ph_uw - scla.ph_scla - iono_corr - tide.ph_tide;
        clear uw scla
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ref_ps=ps_setref;
            
        
    case {'u-o'}
        uw=load(phuwname);
        ph_all = uw.ph_uw;
        % subtract of oscialtor drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);

        clear uw scla
        ref_ps=ps_setref;
        fig_name = 'u-o';
    case {'u-do'}
        uw=load(phuwname);
        scla=load(sclaname);
        ph_all = uw.ph_uw - scla.ph_scla;
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
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
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        
        ref_ps=ps_setref;
        fig_name = 'u-m';
        
     case {'i'}
        iono=load(iononame);
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        fig_name = ['i' fig_name_ica];
        fig_name = [fig_name fig_name_suffix];
        ph_all=iono_corr;
        clear iono iono_corr
        ref_ps=ps_setref;
        
     case {'u-i'}
        uw=load(phuwname);
        iono=load(iononame);
        [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
        fig_name = ['u-i' fig_name_ica];
        fig_name = [fig_name fig_name_suffix];
        ph_all=uw.ph_uw - iono_corr;
        ph_all(:,master_ix)=0;
        clear uw scla
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        ref_ps=ps_setref;
               
    case {'t'}
            tide=load(tidename);
            fig_name = ['t'];
            ph_all= tide.ph_tide;  
            ref_ps=ps_setref;

    case {'a'}
        aps=load(apsname);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['a' fig_name_tca];
        fig_name = [fig_name fig_name_suffix];
        ph_all=aps_corr;
        clear aps aps_corr
        
        if aps_flag~=11
            ref_ps=ps_setref;
        else
            units = ['rad/m^{\alpha}'];
        end
        
   case {'a-o'}
        aps=load(apsname);
        %ph_all=exp(j*ph_scn);
        [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
        fig_name = ['a-o' fig_name_tca];
        fig_name = [fig_name fig_name_suffix];
        ph_all=aps_corr;
        % deramping ifgs
        [ph_all] = ps_deramp(ps,ph_all);
        clear aps aps_corr
        ref_ps=ps_setref;
        
    case {'s'}
        scn=load(scnname);
        %ph_all=exp(j*ph_scn);
        ph_all=scn.ph_scn_slave;
        clear scn
        ref_ps=ps_setref;
        fig_name = 's';
   case {'u-s'}
        uw=load(phuwname);
        scn=load(scnname);
        ph_all=uw.ph_uw;
        ph_all=uw.ph_uw - scn.ph_scn_slave;
        ph_all(:,master_ix)=0;
        clear uw scla
        % subtract of oscilator drift in case of envisat
        ph_all = ph_all-ph_unw_eni_osci;
        
        ref_ps=ps_setref;
        fig_name = 'u-s';
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
        scla=load(sclaname,'K_ps_uw');
        ph_all=scla.K_ps_uw;
        %if exist([pmname,'.mat'],'file')
        %    pm=load(pmname,'K_ps');
        %    ph_all=ph_all+pm.K_ps;
        %end
        clear scla
        ref_ps=ps_setref;
        units='rad/m';
        fig_name = 'd';
    case {'o'}
        if strcmp('n',scla_deramp)
            disp('Warning: scla_deramp flag set to n. Set to y and rerun Step 7 before using the -o plot command.')
            return;
        end
        scla=load(sclaname);
        ph_all=scla.ph_ramp;
        clear scla
        ref_ps=ps_setref;
        fig_name = 'o';
    case {'v'}
        uw=load(phuwname);
        ph_uw=uw.ph_uw;
        % subtract of oscilator drift in case of envisat
        ph_uw = ph_uw-ph_unw_eni_osci;
                
        clear uw
        switch(value_type)
        case {'v'}
            fig_name = 'v';
        case {'v-d'}
            scla=load(sclaname);
            if isfield(scla,'C_ps_uw')
                ph_uw=ph_uw - scla.ph_scla - repmat(scla.C_ps_uw,1,size(ph_uw,2)); % master phase doesn't effect plot, but better for ts plot
            else
                 ph_uw=ph_uw - scla.ph_scla;
            end
            clear scla
            fig_name = 'v-d';
        case {'v-o'}
            ph_uw=ph_uw;
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            clear scla
            fig_name = 'v-o';
            
        case {'v-do'}
            scla=load(sclaname);
            ph_uw=ph_uw - scla.ph_scla;
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            clear scla
            fig_name = 'v-do';

        case {'v-s'}
             scn=load(scnname);
             ph_uw=ph_uw - scn.ph_scn_slave;
             clear scn
             fig_name = 'v-s';
             
        case {'v-so'}
             scn=load(scnname);
             ph_uw=ph_uw - scn.ph_scn_slave;
             clear scn
             % deramping ifgs
             [ph_uw] = ps_deramp(ps,ph_uw);

             fig_name = 'v-so';

         case {'v-dso'}
             scla=load(sclaname);
             scn=load(scnname);
             ph_uw=ph_uw - scla.ph_scla - scn.ph_scn_slave;
             clear scla scn
             % deramping ifgs
             [ph_uw] = ps_deramp(ps,ph_uw);
             fig_name = 'v-dso';

        case {'v-a'}
            aps=load(apsname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['v-a' fig_name_tca];  
            ph_uw=ph_uw - aps_corr;
            clear scla aps aps_corr
            
        case {'v-ao'}
            aps=load(apsname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['v-ao' fig_name_tca];
            ph_uw=ph_uw - aps_corr;
            clear scla aps aps_corr  
            
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
             
        case {'v-da'}
            scla=load(sclaname);
            aps=load(apsname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['v-da' fig_name_tca];
            ph_uw=ph_uw - scla.ph_scla - aps_corr;
            clear scla aps aps_corr
            
        case {'v-dai'}
            scla=load(sclaname);
            aps=load(apsname);
            iono=load(iononame);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
            fig_name = ['v-dai' fig_name_tca ' ' fig_name_ica];
            ph_uw=ph_uw - scla.ph_scla - aps_corr - iono_corr;
            clear scla aps aps_corr iono_corr
         case {'v-dait'}
            scla=load(sclaname);
            aps=load(apsname);
            iono=load(iononame);
            tide=load(tidename);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
            fig_name = ['v-dait' fig_name_tca ' ' fig_name_ica];
            ph_uw=ph_uw - scla.ph_scla - aps_corr -iono_corr - tide.ph_tide;
            clear scla aps aps_corr iono_corr
            
         case {'vt'}
            tide=load(tidename);
            fig_name = ['vt'];
            ph_uw= tide.ph_tide;
            
            
         case {'vi'}
            iono=load(iononame);
            [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
            fig_name = ['vi' fig_name_ica];
            ph_uw= iono_corr; 
            
        case {'v-dao'}
            scla=load(sclaname);
            aps=load(apsname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['v-dao' fig_name_tca];
            ph_uw=ph_uw - scla.ph_scla - aps_corr;
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            
            clear scla aps aps_corr
        case {'v-ds'}
            scla=load(sclaname);
            scn=load(scnname);
            ph_uw=ph_uw - scla.ph_scla - scn.ph_scn_slave;
            clear scla scn
            fig_name = 'v-ds';
        case {'vdrop'}
            fig_name = 'vdrop';
        case {'vdrop-d'}
            scla=load(sclaname);
            ph_uw=ph_uw - scla.ph_scla;
            clear scla
            fig_name = 'vdrop-d';
        case {'vdrop-o'}
            scla=load(sclaname);
            clear scla
            fig_name = 'vdrop-o';
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            
        case {'vdrop-do'}
            scla=load(sclaname);
            ph_uw=ph_uw - scla.ph_scla;
            clear scla
            fig_name = 'vdrop-do';
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
        case {'vdrop-dos'}
             scla=load(sclaname);
             scn=load(scnname);
             ph_uw=ph_uw - scla.ph_scla - scn.ph_scn_slave;
             clear scla scn
             fig_name = 'vdrop-dos';
             % deramping ifgs
             [ph_uw] = ps_deramp(ps,ph_uw);

        case {'vdrop-dao'}
            scla=load(sclaname);
            aps=load(apsname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['vdrop-dao' fig_name_tca];
            ph_uw=ph_uw - scla.ph_scla - aps_corr;
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            
        otherwise
            error('unknown value type')
        end
        
        

        if ts_flag==1 % master AOE doesn't effect v plot, but better for ts plot
            if exist([sclaname '.mat'],'file')==2
                scla=load(sclaname,'C_ps_uw');
                if isfield(scla,'C_ps_uw')
                    ph_uw=ph_uw - repmat(scla.C_ps_uw,1,size(ph_uw,2));
                    clear scla
                else
                   fprintf('Master atmosphere is not subtracted \n') 
                end
            else
               fprintf('Master atmosphere is not subtracted \n') 
            end
        end

        ph_all=zeros(n_ps,1);
        ref_ps=ps_setref;
        
%AH2CHECK 
% % % %         if unwrap_ifg_index(1)~=ps.master_ix & unwrap_ifg_index(end)~=ps.master_ix
% % % %             unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix); % need to include it if not ifgs either side of master
% % % %         end
        %%%% change to always remove the master, not sure why to keep
        %%%% master in for first/last acquisition as the covariance matrix would
        %%%% be having a zeros row and column at master_ix
        unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
        
        if ~isempty(ifg_list)
            unwrap_ifg_index=intersect(unwrap_ifg_index,ifg_list);
            ifg_list=[];
        end
        ph_uw=ph_uw(:,unwrap_ifg_index);
        day=day(unwrap_ifg_index,:);
        
        %ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:),1),n_ps,1);
        ph_uw=ph_uw-repmat(nanmean(ph_uw(ref_ps,:),1),n_ps,1);
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
                m=lscov(G([1:i-1,i+1:end],:),double(ph_uw(:,[1:i-1,i+1:n]))',sm_cov([1:i-1,i+1:end],[1:i-1,i+1:end]));
                ph_all(:,i)=-m(2,:)'*365.25/4/pi*lambda*1000; 
            end
        else     
            m=lscov(G,double(ph_uw'),sm_cov);
            ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
        end

        try
            save mean_v m
        catch
            save ~/mean_v m
            fprintf('Warning: Read access only, velocities saved in home directory\n')
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
        % subtract of oscialtor drift in case of envisat
        ph_uw = ph_uw-ph_unw_eni_osci;
        
        clear phuw
        switch(value_type)
        case('v')
            fig_name = 'v';
        case('v-d')
            scla=load(sclasbname);
            ph_uw=ph_uw - scla.ph_scla;
            clear scla
            fig_name = 'v-d';
               
        case {'v-a'}
            aps=load(apssbname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['v-a' fig_name_tca];
            ph_uw=ph_uw - aps_corr;
            clear scla aps aps_corr
            
         case {'v-ao'}
            aps=load(apssbname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['v-ao' fig_name_tca];
            ph_uw=ph_uw - aps_corr;
            clear scla aps aps_corr
            
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            
        case {'v-da'}
            scla=load(sclasbname);
            aps=load(apssbname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['v-da' fig_name_tca];
            ph_uw=ph_uw - scla.ph_scla - aps_corr;
            clear scla aps aps_corr
        case {'v-dai'}
            scla=load(sclasbname);
            aps=load(apssbname);
            iono=load(ionosbname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
            fig_name = ['v-dai' fig_name_tca ' ' fig_name_ica];
            ph_uw=ph_uw - scla.ph_scla - aps_corr -iono_corr;
            clear scla aps aps_corr iono_corr
                  
            
        case {'v-dat'}
            scla=load(sclasbname);
            aps=load(apssbname);
            tide=load(tidesbname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['v-dat' fig_name_tca];
            ph_uw=ph_uw - scla.ph_scla - aps_corr- tide.ph_tide;
            clear scla aps aps_corr  
            
          case {'v-dait'}
            scla=load(sclasbname);
            aps=load(apssbname);
            iono=load(ionosbname);
            tide=load(tidesbname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            [iono_corr,fig_name_ica] = ps_plot_ica(iono,iono_flag);
            fig_name = ['v-dait' fig_name_tca ' ' fig_name_ica];
            ph_uw=ph_uw - scla.ph_scla - aps_corr -iono_corr - tide.ph_tide;
            clear scla aps aps_corr iono_corr
              
            
        case {'v-dao'}
            scla=load(sclasbname);
            aps=load(apssbname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['v-dao' fig_name_tca];
            ph_uw=ph_uw - scla.ph_scla - aps_corr;
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            
            clear scla aps aps_corr
        case('v-o')
            ph_uw=ph_uw;
            fig_name = 'v-o';
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            
        case('v-do')
            scla=load(sclasbname);            
            ph_uw=ph_uw - scla.ph_scla ;
            ph_uw = ps_deramp(ps,ph_uw);
            clear scla
            fig_name = 'v-do';
        case {'vdrop'}
            fig_name = 'vdrop';
        case {'vdrop-d'}
            scla=load(sclasbname);
            ph_uw=ph_uw - scla.ph_scla;
            clear scla
            fig_name = 'vdrop-d';
        case {'vdrop-o'}
            scla=load(sclasbname);
            ph_uw=ph_uw;
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            clear scla
            fig_name = 'vdrop-o';
        case {'vdrop-do'}
            scla=load(sclasbname);
            ph_uw=ph_uw - scla.ph_scla;
            clear scla
            fig_name = 'vdrop-do';
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            
        case {'vdrop-dao'} 
            scla=load(sclasbname);
            aps=load(apssbname);
            [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag);
            fig_name = ['vdrop-dao' fig_name_tca];
            ph_uw=ph_uw - scla.ph_scla - aps_corr;
            % deramping ifgs
            [ph_uw] = ps_deramp(ps,ph_uw);
            
            clear scla aps aps_corr
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
        % fix for aps data that has a nan value in the reference area
        if sum(sum(isnan(ph_uw)))>0
            for kk=1:size(ph_uw,2)
                    ph_uw(:,kk)=ph_uw(:,kk)-repmat(mean(ph_uw(~isnan(ph_uw(:,kk)),kk),1),n_ps,1);
            end
        else
            ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:),1),n_ps,1);
        end

        G=[ones(size(ifgday_ix(:,1))),day(ifgday_ix(:,2))-day(ifgday_ix(:,1))];
        lambda=getparm('lambda');

        if length(value_type)>4 & strcmpi(value_type(1:5),'vdrop') 
            ph_all=zeros(size(ph_uw));
            n=size(ph_uw,2);
            for i=1:n
                m=lscov(G([1:i-1,i+1:end],:),double(ph_uw(:,[1:i-1,i+1:n])'),sb_cov([1:i-1,i+1:end],[1:i-1,i+1:end]));
                ph_all(:,i)=-m(2,:)'*365.25/4/pi*lambda*1000; 
            end
        else
            m=lscov(G,double(ph_uw'),sb_cov);
            ph_all=-m(2,:)'*365.25/4/pi*lambda*1000; 
        end
        try
            save mean_v m
        catch
            save ~/mean_v m
            fprintf('Warning: Read access only, velocities saved in home directory\n')
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
%         fig_name_tca=ps_mean_v_cov(ifg_list,800,value_type(4:end),use_small_baselines,aps_flag,reference_flag);
        fig_name_tca=ps_mean_v(ifg_list,1500,value_type(4:end),use_small_baselines,aps_flag);

        ifg_list=[];
        mv=load(meanvname);
        ph_all=mv.mean_v_std;
        units='mm/yr';
        fig_name = [value_type fig_name_tca];
    otherwise
        error('unknown value type')
    end
else
	ph_all = value_type;
    ref_ps=ps_setref;
	fig_name = 'data';
    value_type = 'data';
    if size(ph_all,2)~=ps.n_ifg
        ifg_list=[];
    end
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
            if isnan(mean_ph(i))
                mean_ph(i)=0;
                fprintf(['Interferogram (' num2str(i) ') does not have a reference area\n'])
            end
        end
        clear i
        ph_disp=ph_disp-repmat(mean_ph,n_ps,1);
    end
    phsort=sort(ph_disp(~isnan(ph_disp)));
    
    if isempty(lims)
        maxph=phsort(round(length(phsort)*.999));
        minph=phsort(ceil(length(phsort)*.001));
        lims=[minph,maxph];
        if isempty(lims)==1
            fprintf(['Interferograms do not contain data.\n'])
        end
        
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
        clear i
        ph_disp=ph_disp.*conj(repmat(mean_ph,n_ps,1));
    end
    lims=[-pi,pi];

end

if plot_flag==-1
    %savename=['~/ps_plot_',value_type]
    h_fig=[];
    lims=[];
    ifg_data_RMSE=[];
    h_axes_all=[];
    savename=['ps_plot_',value_type]
    try
       stamps_save(savename,ph_disp,ifg_list)
    catch
       stamps_save(['~/',savename],ph_disp,ifg_list)
       fprintf('Warning: Read access only, values in home directory instead\n')
    end
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
  ifg_data_RMSE=NaN([size(ph_disp,2) 1]) ;
  if size(ifg_list,1)>1
    ifg_list = ifg_list';
  end
  for i=ifg_list
    i_im=i_im+1;
    if n_ifg_plot>1
        h_axes = axes('position',[imX(i_im),imY(i_im),h_x,h_y]);
        set(h_axes,'Xcolor',[1 1 1]);			% [DB] remove the contours of the axes figure
        set(h_axes,'Ycolor',[1 1 1]);
        
        h_axes_all(i)=h_axes;
        clear h_axes
    else
        h_axes_all(i)=gca;
    end
    % check if external data is requested to be plotted
    if ext_data_flag==1
        if exist(ext_data_path,'file')==2
            % this is an individual file identical for each ifgs
            loadname=ext_data_path;
        else      
            if size(day,1)==size(ph_all,2) && strcmpi(small_baseline_flag,'y')
                % This is for external data for a SB to a SM inverted network
                if aps_band_flag==1
                    fprintf('Still to do')
                    keyboard 
                else                   
                    % checking if there is data for this interferogram
                    loadname = [ext_data_path filesep  datestr(ps.day(i,1),'yyyymmdd') '_' datestr(ps.master_day,'yyyymmdd') '.mat'];
                end
 
            else
                % This is external data for a SM or SB network. It is not a
                % SB to SM inverted network.
                % The path is specified
                if aps_band_flag==1
                    % checking if there is data for this interferogram
                    loadname = [ext_data_path filesep datestr(ps.ifgday(ifg_number,1),'yyyymmdd') '_' datestr(ps.ifgday(ifg_number,2),'yyyymmdd') '.mat'];
                else
                    % checking if there is data for this interferogram
                    loadname = [ext_data_path filesep datestr(ps.ifgday(i,1),'yyyymmdd') '_' datestr(ps.ifgday(i,2),'yyyymmdd') '.mat'];
                end
 
            end
            
            
            
           
        end
        if exist(loadname,'file')==2
             % loading of the data 
             ext_data = load(loadname);
             
             
             % Checking if there is a ph_disp variable stored in the files 
             if isfield(ext_data,{'ph_disp'})
                 % Correct the data with respect to the minimization of the
                 % residual between this dataset and the data
                 % convert to a local reference
                 n_data_points = size(ext_data.lonlat,1);
                 data_xy=llh2local(ext_data.lonlat',ps.ll0)*1000;
                 ps_xy=llh2local(ps.lonlat',ps.ll0)*1000;

                 dist_sq=(repmat(ps_xy(1,:)',1,n_data_points) -repmat(data_xy(1,:),ps.n_ps,1)).^2+(repmat(ps_xy(2,:)',1,n_data_points) -repmat(data_xy(2,:),ps.n_ps,1)).^2; 
                 for kk=1:n_data_points      
                     ref_data= (dist_sq(:,kk)<=ref_radius_data^2);
                     if sum(ref_data)>0
                        ifg_mean_data_point(kk,1) = mean(ph_disp(ref_data,i_im));
                     else
                        ifg_mean_data_point(kk,1) = NaN;  
                     end
                 end
                 clear ref_data
                 if ~isempty(ifg_mean_data_point) && sum(isnan(ifg_mean_data_point))~=n_data_points
                     % mean residual 
                     ext_data.mean_residual = mean(ext_data.ph_disp(~isnan(ifg_mean_data_point),1)-ifg_mean_data_point(~isnan(ifg_mean_data_point)));
                     % correct data phases
                     ext_data.ph_disp(:,1) = ext_data.ph_disp(:,1)-ext_data.mean_residual;
                     % computing the RMSE
                     ext_data.residual = ext_data.ph_disp(~isnan(ifg_mean_data_point),1)-ifg_mean_data_point(~isnan(ifg_mean_data_point));
                     ext_data.RMSE = sqrt(mean((ext_data.residual).^2));
                     clear ifg_mean_data_point
                     ifg_data_RMSE(i)=ext_data.RMSE;
                 else
                    fprintf(['No observation within the ref_radius_data, increase the size. \n']) 
                    fprintf(['External data not plotted for this interferogram. \n']) 
                    ext_data = [];
                 end
                 clear mean_data_point_residual

    %              % correct the data with respect to the reference area
    %              point_ref = ps_setref(ext_data);
    %              ext_data.ph_disp(:,1) = ext_data.ph_disp(:,1)-mean(ext_data.ph_disp(point_ref,1));
            end

        else
           ext_data = []; 
        end
    else
       ext_data = []; 
    end
    ps_plot_ifg(ph_disp(:,i_im),plot_flag,lims,lon_rg,lat_rg,ext_data);
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
        if textsize~=0 & size(day,1)==size(ph_all,2) & aps_band_flag==0 && (strcmpi(small_baseline_flag,'n') | strcmpi(value_type(1),'s'))
            % text for SM ifgs
            t=text(x_t,y_t,[datestr(day(i),'dd mmm yyyy')]);
            set(t,'fontweight','bold','color',textcolor,'fontsize',abs(textsize))
        elseif textsize~=0 & size(day,1)==size(ph_all,2) & aps_band_flag==0 && forced_sm_flag==1
            % text for SM ifgs when inverting from SBAS
            t=text(x_t,y_t,[datestr(day(i),'dd mmm yyyy')]);
            set(t,'fontweight','bold','color',textcolor,'fontsize',abs(textsize))
        elseif textsize~=0 & ps.n_ifg==size(ph_all,2) & aps_band_flag==0 && strcmpi(small_baseline_flag,'y')
            % text for SB ifgs
            t=text(x_t,y_t,['      ifg ' num2str(i)]);
            set(t,'fontweight','bold','color',textcolor,'fontsize',abs(textsize))
        elseif textsize~=0 & size(bands,1)==size(ph_all,2) & aps_band_flag==1
            % text for band filtered data
            t=text(x_t,y_t,[num2str(round((bands(i,1)/100))*100/1000) ' - ' num2str(round((bands(i,2)/100))*100/1000) ' km']);
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
                set(h,'xtick',[xlim(2)-64,xlim(2)],'Xticklabel',plotlims,'xcolor','k','ycolor',textcolor2,'fontweight','bold','color',textcolor2,'FontSize',abs(textsize))
                h=xlabel(h,units);
                pos=get(h,'position');
                pos(2)=pos(2)/2.2;
                set(h,'position',pos,'FontSize',abs(textsize));
        else
                set(h,'ytick',[ylim(2)-64,ylim(2)],'yticklabel',plotlims,'xcolor','k','ycolor',textcolor2,'fontweight','bold','color',textcolor2,'FontSize',abs(textsize))
                set(get(h,'ylabel'),'String',units,'FontSize',abs(textsize))  

        end
    end
  end
end

if ext_data_flag==1
    for k=1:size(ifg_data_RMSE,1)
        if k==1
            fprintf(['\nRMSE between interferogram(s) and external data \n'])
        end
        
        if aps_band_flag==1
            fprintf(['Band' num2str(k) ' : ' num2str(ifg_data_RMSE(k)) '\n']);
        else
            if size(ifg_data_RMSE,1)==ps.n_ifg
                fprintf(['ifg ' num2str(k) ' \t ' datestr(ps.ifgday(k,1),'yyyymmdd') '-' datestr(ps.ifgday(k,2),'yyyymmdd') ' \t ' num2str(ifg_data_RMSE(k)) '\n']);
            else
                 fprintf(['dataset ' num2str(k) ' \t ' num2str(ifg_data_RMSE(k)) '\n']);               
            end
        end
    end
end


if exist('plot_color_scheme_old','var')==1
   setparm('plot_color_scheme',plot_color_scheme_old) 
end

fprintf('Color Range: %g to %g %s\n',lims,units)
if ts_flag == 1
  figure(h_fig);
  clear all % clean up to save memory
  
  % Place button for new TS plots
  fPosition=get(gcf,'position');
  %  pos_new = [pos(1) pos(2)+0.1 pos(3) pos(4)-0.1]
  
  % new TS plot button
  mButton=uicontrol('Style', 'pushbutton', 'Callback', 'clear ph_uw; ts_plot',...
    'String','TS plot', 'Position', [150 30 90 20] ); % if exist don't create
  %                                      left bottom width height
  %mButtonposition=get(mButton,'Position')
 
  % new TS plot button for selecting two points and difference
  mButton=uicontrol('Style', 'pushbutton', 'Callback', 'clear ph_uw; ts_plotdiff',...
    'String','TS double diff.', 'Position', [470 30 90 20] ); % 
  
  % Radius Text boxes
  mTextBox=uicontrol('Style', 'edit','String','radius (m) ', 'Position',...
      [320 30 90 20] );

  mEditBox=uicontrol('Style', 'edit','String','100', 'Position',...
      [410 30 30 20],'BackgroundColor',[1 1 1] );
%  ts_plot   % select a point then plot, for the first time.
end




    
