function [ph_uw,msd]=uw_3d(ph,xy,day,ifgday_ix,bperp,options)
%UW_3D unwrap phase time series (single or multiple master)
%   PH_UW = UW_3D(PH,XY,DAY,IFGDAY_IX,BPERP,OPTIONS)
%
%   PH  = N x M matrix of wrapped phase values (real phase or complex phasor)
%        where N is number of pixels and M is number of interferograms
%   XY  = N x 2 matrix of coordinates in metres
%        (optional extra column, in which case first column is ignored)
%   DAY = vector of image acquisition dates in days, relative to master
%   IFGDAY_IX = M x 2 matrix giving index to master and slave date in DAY
%        for each interferogram 
%   BPERP  = M x 1 vector giving perpendicular baselines 
%   OPTIONS structure optionally containing following fields:
%      la_flag    = look angle error estimation flag (def='y')
%      scf_flag   = spatial cost function estimation flag (def='y')
%      master_day = serial date number of master (used for info msg only, def=0)
%      grid_size  = size of grid in m to resample data to (def=5)
%      prefilt_win = size of prefilter window in resampled grid cells (def=16)
%      time_win   = size of time filter window in days (def=365)
%      unwrap_method (def='3D' for single master, '3D_FULL' otherwise)
%      goldfilt_flag = Goldstein filtering, 'y' or 'n' (def='n')
%      gold_alpha    = alpha value for Goldstein filter (def=0.8)
%      lowfilt_flag  = Low pass filtering, 'y' or 'n' (def='n')
%      n_trial_wraps = no. phase cycles poss between neighbouring points due to la error (def=6)
%      temp          = optional M x 1 vector of temperature difference for each ifg (def=[])
%      n_temp_wraps  = no. phase cycles poss between neighbouring points due to temp diff (def=2)
%      variance      = N x 1 matrix of variance values           
%   PH_UW = unwrapped phase
%
%   Andy Hooper, Jun 2007

% ============================================================
% 08/2009 AH: Goldstein alpha value added to options
% 01/2012 AH: Changes for easier running stand-alone
% 01/2012 AH: Add bperp for new look angle error estimation
% 02/2012 AH: Add new method 3D_NEW
% 11/2012 AH: Several new options
% 02/2014 AH: Add predefined ph_uw option
% ============================================================
tic;
if nargin<4
    help uw_3d
    error('not enough arguments')
end

if nargin<5
    bperp=[];
end

if nargin<6
    options=struct;
end

%if isempty(ifgday_ix) | length(unique(ifgday_ix(:,1)))==1
if isempty(ifgday_ix) 
    single_master_flag=1;
else
    single_master_flag=0;
end

valid_options={'la_flag','scf_flag','master_day','grid_size','prefilt_win','time_win','unwrap_method','goldfilt_flag','lowfilt_flag','gold_alpha','n_trial_wraps','temp','n_temp_wraps','max_bperp_for_temp_est','variance','ph_uw_predef'};
invalid_options=setdiff(fieldnames(options),valid_options);
if length(invalid_options)>0
    error(['"',invalid_options{1}, '" is an invalid option'])
end

if ~isfield(options,'master_day')
    options.master_day=0;
end

if ~isfield(options,'grid_size')
    options.grid_size=5;
end

if ~isfield(options,'prefilt_win')
    options.prefilt_win=16;
end

if ~isfield(options,'time_win')
    options.time_win=365;
end

if ~isfield(options,'unwrap_method')
    if single_master_flag==1
        options.unwrap_method='3D';
    else
        options.unwrap_method='3D_FULL';
    end
end

if ~isfield(options,'goldfilt_flag')
    options.goldfilt_flag='n';
end

if ~isfield(options,'lowfilt_flag')
    options.lowfilt_flag='n';
end

if ~isfield(options,'gold_alpha')
    options.gold_alpha=0.8;
end

if ~isfield(options,'n_trial_wraps')
    options.n_trial_wraps=6;
end

if ~isfield(options,'la_flag')
    options.la_flag='y';
end

if ~isfield(options,'scf_flag')
    options.scf_flag='y';
end

if ~isfield(options,'temp')
    options.temp=[];
else
    if ~isempty(options.temp) & length(options.temp)~=size(ph,2)
        error('options.temp must be M x 1 vector where M is no. of ifgs')
    end
end

if ~isfield(options,'n_temp_wraps')
    options.n_temp_wraps=2;
end

if ~isfield(options,'max_bperp_for_temp_est')
    options.max_bperp_for_temp_est=100;
end

if ~isfield(options,'variance')
    options.variance=[];
end

if ~isfield(options,'ph_uw_predef')
    options.ph_uw_predef=[];
end

if size(xy,2)==2
   xy(:,2:3)=xy(:,1:2);
end

if size(day,1)==1
    day=day';
end
    
if strcmpi(options.unwrap_method,'3D') | strcmpi(options.unwrap_method,'3D_NEW') 
    if length(unique(ifgday_ix(:,1)))==1
        options.unwrap_method='3D_FULL';
    else
        options.lowfilt_flag='y';
    end
end

uw_grid_wrapped(ph,xy,options.grid_size,options.prefilt_win,options.goldfilt_flag,options.lowfilt_flag,options.gold_alpha,options.ph_uw_predef);
clear ph
uw_interp;
%if single_master_flag==1
%    uw_unwrap_space_time(day,options.unwrap_method,options.time_win,options.master_day,options.la_flag,bperp,options.n_trial_wraps,options.prefilt_win,options.scf_flag,options.temp,options.n_temp_wraps);
%else
    uw_sb_unwrap_space_time(day,ifgday_ix,options.unwrap_method,options.time_win,options.la_flag,bperp,options.n_trial_wraps,options.prefilt_win,options.scf_flag,options.temp,options.n_temp_wraps,options.max_bperp_for_temp_est);
%end
uw_stat_costs(options.unwrap_method,options.variance);
[ph_uw,msd]=uw_unwrap_from_grid(xy,options.grid_size);

