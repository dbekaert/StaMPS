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
%        for each interferogram (can be empty for single master time series)
%   BPERP  = M x 1 vector giving perpendicular baselines 
%   OPTIONS structure optionally containing following fields:
%      la_flag    = look angle error estimation flag (def='y')
%      master_day = serial date number of master (used for info msg only, def=0)
%      grid_size  = size of grid in m to resample data to (def=5)
%      prefilt_win = size of prefilter window in resampled grid cells (def=16)
%      time_win   = size of time filter window in days (def=90)
%      unwrap_method (def='3D' for single master, '3D_QUICK' otherwise)
%      goldfilt_flag = Goldstein filtering, 'y' or 'n' (def='n')
%      gold_alpha    = alpha value for Goldstein filter (def=0.8)
%      lowfilt_flag  = Low pass filtering, 'y' or 'n' (def='n')
%      n_trial_wraps = for la_flag='y', num phase cycles poss over bperp range (def=6)
%   PH_UW = unwrapped phase
%
%   Andy Hooper, Jun 2007

% ============================================================
% 08/2009 AH: Goldstein alpha value added to options
% 01/2012 AH: Changes for easier running stand-alone
% 01/2012 AH: Add bperp for new look angle error estimation
% 02/2012 AH: Add new method 3D_NEW
% ============================================================

if nargin<3
    help uw_3d
    error('not enough arguments')
end

if nargin<4
    ifgday_ix=[];
end

if nargin<5
    bperp=[];
end

if nargin<6
    options=struct;
end

if isempty(ifgday_ix) | length(unique(ifgday_ix(:,1)))==1
    single_master_flag=1;
else
    single_master_flag=0;
end

valid_options={'la_flag','master_day','grid_size','prefilt_win','time_win','unwrap_method','goldfilt_flag','lowfilt_flag','gold_alpha'};
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
    options.time_win=90;
end

if ~isfield(options,'unwrap_method')
    if single_master_flag==1
        options.unwrap_method='3D_NEW';
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

if size(xy,2)==2
   xy(:,2:3)=xy(:,1:2);
end

if size(day,1)==1
    day=day';
end

if strcmpi('unwrap_method','3D') & single_master_flag~=1
    options.lowfilt_flag='y';
end

uw_grid_wrapped(ph,xy,options.grid_size,options.prefilt_win,options.goldfilt_flag,options.lowfilt_flag,options.gold_alpha);
uw_interp;
if single_master_flag==1
    uw_unwrap_space_time(day,options.unwrap_method,options.time_win,options.master_day,options.la_flag,bperp,options.n_trial_wraps,options.prefilt_win);
else
    uw_sb_unwrap_space_time(day,ifgday_ix,options.unwrap_method,options.time_win,options.la_flag,bperp,options.n_trial_wraps,options.prefilt_win);
end
uw_stat_costs(options.unwrap_method);
[ph_uw,msd]=uw_unwrap_from_grid(ph,xy,options.grid_size);

