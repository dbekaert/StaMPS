function ps_parms_default()
%PS_PARMS_DEFAULT set parms to default value if not already set
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   04/2008 AH: SB default processing corrected
%   06/2009 AH: scla_deramp added 
%   08/2009 AH: unwrap_gold_alpha added
%   09/2009 AH: merge_resample_size added 
%   09/2009 AH: select_method added
%   09/2009 AH: density_rand added
%   11/2009 AH: pixel_aspect_ratio removed
%   01/2010 AH: add check for write permission before saving
%   02/2010 AH: weed_alpha replaced by weed_time_win
%   02/2010 AH: weed_max_noise added
%   02/2010 AH: plot_pixel_size replaced by plot_scatterer_size (in m)
%   02/2010 AH: plot_dem_posting added
%   02/2010 AH: ref_centre_lonlat and ref_radius added
%   02/2010 AH: unwrap_ifg_index replaced by drop_ifg_index
%   02/2010 AH: weed_time_win and unwrap_time_win default changed to 730
%   03/2010 AH: add logging
%   09/2010 MA: oversampling factor
%   09/2010 AH: plot_pixels_scatterer added
%   11/2010 AH: (sb_)recalc_index replaced by (sb_)scla_drop_index
%   06/2011 AH: add weed_neighbours 
%   12/2012 AH: add insar_processor
%   06/2013 DB: add parameters to include tropopshere in scla estimation
%   01/2014 DB: Add wavelength and heading to parameters list. This in case
%               not all stamps steps are ran.
%   01/2014 DB: Add a multi-core option for step 1-5
%   02/2014 AH: add unwrap_hold_good_values
%   08/2014 DB: explicit search if there is not a processor.txt file 
%   05/2015 AH: add gamma_max_iterations 
%   09/2015 AH: change default weed_neighbours = 'n' 
%   09/2015 AH: change default max_topo_err = 20 
%   10/2016 DB: Change other processor type to a warning not an error.
%   ======================================================================


parmfile='parms.mat';
parent_flag=0;

if exist('./parms.mat','file')
    parms=load(parmfile);
elseif exist('../parms.mat','file')
    parmfile='../parms.mat';
    parms=load(parmfile);
    parent_flag=1;
else
    parms=struct('Created',date);
    parms.small_baseline_flag='n'; 
end

parmfields_before=fieldnames(parms);
num_fields=size(parmfields_before,1);

if ~isfield(parms,'max_topo_err')
    parms.max_topo_err=20;
end

if ~isfield(parms,'quick_est_gamma_flag')
    parms.quick_est_gamma_flag='y';
end

if ~isfield(parms,'select_reest_gamma_flag')
    parms.select_reest_gamma_flag='y';
end

if ~isfield(parms,'filter_grid_size')
    parms.filter_grid_size=50;
end

if ~isfield(parms,'filter_weighting')
    parms.filter_weighting='P-square'; % filter weighting strategy
end

if ~isfield(parms,'gamma_change_convergence')
    parms.gamma_change_convergence=0.005; % change in change in gamma that signals convergence
end

if ~isfield(parms,'gamma_max_iterations')
    parms.gamma_max_iterations=3; % maximum number of iterations for gamma estimation 
end

if ~isfield(parms,'slc_osf')
    parms.slc_osf=1;  % [MA] SLC oversampling factor
end

if ~isfield(parms,'clap_win')
    parms.clap_win=32;
end

if ~isfield(parms,'clap_low_pass_wavelength')
    parms.clap_low_pass_wavelength=800;
end

if ~isfield(parms,'clap_alpha')
    parms.clap_alpha=1;
end

if ~isfield(parms,'clap_beta')
    parms.clap_beta=0.3;
end

if ~isfield(parms,'select_method')
    parms.select_method='DENSITY'; %'DENSITY' or 'PERCENT'
end
 
if ~isfield(parms,'density_rand')
    if strcmpi(parms.small_baseline_flag,'y')
        parms.density_rand=2;  % Random-phase pixels per km2 (before weeding)
    else
        parms.density_rand=20; % Random-phase pixels per km2 (before weeding) 
    end
end

if ~isfield(parms,'percent_rand')
    if strcmpi(parms.small_baseline_flag,'y')
        parms.percent_rand=1; % Percent random-phase pixels (before weeding) 
    else
        parms.percent_rand=20;% Percent random-phase pixels (before weeding) 
    end
end

if ~isfield(parms,'gamma_stdev_reject')
    parms.gamma_stdev_reject=0;
end

if isfield(parms,'weed_alpha')
    parms=rmfield(parms,'weed_alpha');
end

if ~isfield(parms,'weed_time_win')
    parms.weed_time_win=730;    % weeding smoothing window alpha (days)
end

if ~isfield(parms,'weed_max_noise')
    parms.weed_max_noise=inf;   % maximum arc noise in any ifg (rad)
end

if ~isfield(parms,'weed_standard_dev')
    if strcmpi(parms.small_baseline_flag,'y')
        parms.weed_standard_dev=inf;
    else 
        parms.weed_standard_dev=1.0;
    end
end

if ~isfield(parms,'weed_zero_elevation')
    parms.weed_zero_elevation='n';
end

if ~isfield(parms,'weed_neighbours')
    parms.weed_neighbours='n';
end

if ~isfield(parms,'unwrap_method')
    if strcmpi(parms.small_baseline_flag,'y')
        parms.unwrap_method='3D_QUICK';  
    else
        parms.unwrap_method='3D';
    end
end

if ~isfield(parms,'unwrap_patch_phase')
    parms.unwrap_patch_phase='n';
end

if isfield(parms,'unwrap_ifg_index')
    try
       ps=load('ps2.mat');
    catch
       try
          ps=load('ps1.mat');
       catch
       end
    end
    if exist('ps','var') & ~strcmpi(parms.unwrap_ifg_index,'all')
       parms.drop_ifg_index=setdiff([1:ps.n_ifg],parms.unwrap_ifg_index);  
    end
    parms=rmfield(parms,'unwrap_ifg_index');
    num_fields=0;
end

if ~isfield(parms,'drop_ifg_index')
    parms.drop_ifg_index=[];  
end

if ~isfield(parms,'unwrap_la_error_flag')
    parms.unwrap_la_error_flag='y';
end

if ~isfield(parms,'unwrap_spatial_cost_func_flag')
    parms.unwrap_spatial_cost_func_flag='n';
end

if ~isfield(parms,'unwrap_prefilter_flag')
    parms.unwrap_prefilter_flag='y';
end

if ~isfield(parms,'unwrap_grid_size')
    parms.unwrap_grid_size=200; % prefilter grid size
end

if ~isfield(parms,'unwrap_gold_n_win')
    parms.unwrap_gold_n_win=32; % prefilter goldstein filtering window size
end

if ~isfield(parms,'unwrap_alpha')
    parms.unwrap_alpha=8;    % unwrapping smoothing window alpha
end

if ~isfield(parms,'unwrap_time_win')
    parms.unwrap_time_win=730;    % unwrapping smoothing window alpha
end

if ~isfield(parms,'unwrap_gold_alpha')
    parms.unwrap_gold_alpha=0.8;  % unwrapping goldstein filter alpha
end

if ~isfield(parms,'unwrap_hold_good_values')
    parms.unwrap_hold_good_values='n';       % fix unwrapped values judged to be good
end

if isfield(parms,'recalc_index')
    try
       ps=load('ps2.mat');
    catch
       try
          ps=load('ps1.mat');
       catch
       end
    end
    if exist('ps','var') & ~strcmpi(parms.recalc_index,'all')
        if strcmpi(parms.small_baseline_flag,'y')
           parms.scla_drop_index=setdiff([1:ps.n_image],parms.recalc_index);  
        else
           parms.scla_drop_index=setdiff([1:ps.n_ifg],parms.recalc_index);  
        end
    end
    parms=rmfield(parms,'recalc_index');
    if isfield(parms,'sb_recalc_index')
        if exist('ps','var') & ~strcmpi(parms.sb_recalc_index,'all')
            parms.sb_scla_drop_index=setdiff([1:ps.n_ifg],parms.sb_recalc_index);  
        end
        parms=rmfield(parms,'sb_recalc_index');
    end
    num_fields=0;
end

if ~isfield(parms,'scla_drop_index')
    parms.scla_drop_index=[];  
end

if ~isfield(parms,'scn_wavelength')
    parms.scn_wavelength=100; % spatially correlated noise wavelength
end

if ~isfield(parms,'scn_time_win')
    parms.scn_time_win=365; % 1 year time window
end

if ~isfield(parms,'scn_deramp_ifg')
    parms.scn_deramp_ifg=[];  % deramp these ifgs and add to estimate of scn
end

if ~isfield(parms,'scn_kriging_flag')
    parms.scn_kriging_flag='n';  % use kriging, 'y' or 'n'
end

if ~isfield(parms,'ref_lon')
    parms.ref_lon=[-inf,inf];  % low and high longitude for ref ps
end

if ~isfield(parms,'ref_lat')
    parms.ref_lat=[-inf,inf];  % low and high latitude for ref ps
end

if ~isfield(parms,'ref_centre_lonlat')
    parms.ref_centre_lonlat=[0,0];  % centre lon/lat for ref ps
end

if ~isfield(parms,'ref_radius')
    parms.ref_radius=inf;  % radius from centre for ref ps
end

if ~isfield(parms,'ref_velocity')
    parms.ref_velocity=0;  % velocity of reference point if known
end

if ~isfield(parms,'n_cores')
    parms.n_cores=1;  % n_cores for the muti-core option (step 1-5)
end


if ~isfield(parms,'plot_dem_posting')
    parms.plot_dem_posting=90; 
end

if isfield(parms,'plot_pixel_size')
    parms.plot_scatterer_size=parms.plot_pixel_size*25; 
    num_fields=0;
    parms=rmfield(parms,'plot_pixel_size');
end

if ~isfield(parms,'plot_scatterer_size')
    parms.plot_scatterer_size=120; 
end

if ~isfield(parms,'plot_pixels_scatterer')
    parms.plot_pixels_scatterer=3; 
end

if ~isfield(parms,'plot_color_scheme')
    parms.plot_color_scheme='inflation'; 
end

if isfield(parms,'pixel_aspect_ratio')
    parms=rmfield(parms,'pixel_aspect_ratio');
end

if ~isfield(parms,'shade_rel_angle')
    parms.shade_rel_angle=[90,45]; % look angle for dem shaded relief
end

if ~isfield(parms,'lonlat_offset')
    parms.lonlat_offset=[0,0]; % offset of PS in degrees from dem
end

if ~isfield(parms,'merge_resample_size')
    if strcmpi(parms.small_baseline_flag,'y')
        parms.merge_resample_size=100; % grid size (in m) to resample to during merge of patches
    else
        parms.merge_resample_size=0;   % (0=no resampling)
    end
end

if ~isfield(parms,'merge_standard_dev')
    parms.merge_standard_dev=inf;
end

if ~isfield(parms,'scla_method')
    parms.scla_method='L2'; % method for estmating SCLA, L1- or L2-norm
end

if ~isfield(parms,'scla_deramp')
    parms.scla_deramp='n'; % estimate an orbital ramp before SCLA
end

lambdaname=['lambda.1.in'];         % wavelength
if ~isfield(parms,'lambda')
    if ~exist(lambdaname,'file')
        lambdaname= ['../',lambdaname];
    end
    if ~exist(lambdaname,'file')
        lambdaname= ['../',lambdaname];
    end
    if ~exist(lambdaname,'file')
        parms.lambda=NaN; % Add wavelength
    else
        lambda=load(lambdaname);
        parms.lambda=lambda; % Add wavelength
    end
end

headingname=['heading.1.in'];       % satellite heading
if ~isfield(parms,'heading')
    if ~exist(headingname,'file')
        headingname= ['../',headingname];
    end
    if ~exist(headingname,'file')
        headingname= ['../',headingname];
    end
    if ~exist(headingname,'file')
        parms.heading=NaN; % Add heading
    else
        heading=load(headingname);
        parms.heading=heading; % Add heading
    end
end

if ~isfield(parms,'scla_deramp')
    parms.scla_deramp='n'; % estimate an orbital ramp before SCLA
end

if ~isfield(parms,'sb_scla_drop_index')
    if strcmpi(parms.small_baseline_flag,'y')
        parms.sb_scla_drop_index=[];
    end
end

if ~isfield(parms,'insar_processor')
    processor_file = 'processor.txt';
    if exist(processor_file,'file')~=2
       if exist(['..' filesep processor_file],'file')==2
           processor_file = ['..' filesep processor_file];
           if exist(['..' filesep processor_file],'file')==2
               processor_file = ['..' filesep processor_file];
           end
       end
    end 
           
   if exist(processor_file,'file')~=2
       parms.insar_processor='doris'; % 
   else
       processor = fileread(processor_file);
       processor = strtrim(processor);
       parms.insar_processor=processor; % 

       if ~strcmpi(processor,'gamma') & ~strcmpi(processor,'doris')
           fprintf('WARNING: This processor is not supported (doris and gamma)')
       end
   end
end

if ~isfield(parms,'subtr_tropo')
    parms.subtr_tropo='n';         % remove tropospheric estimate
end

if ~isfield(parms,'tropo_method')
    parms.tropo_method='a_l';       % method for tropopsheric estimate
end


parmfields=fieldnames(parms);
if size(parmfields,1)~=num_fields
    try
        save(parmfile,'-struct','parms')
        for i=1:size(parmfields,1)
            if isempty(strmatch(parmfields{i},parmfields_before))
               parmname=parmfields{i};
               value=getfield(parms,parmname);
               if isempty(value)
                   value='[]';
               end
               if isnumeric(value)
                   logit([parmname,' = ',num2str(value)],0,parent_flag);
               else
                   logit([parmname,' = ',value],0,parent_flag);
               end
            end
        end

    catch
        fprintf('Warning: missing parameters could not be updated (no write access)\n')
    end
end
