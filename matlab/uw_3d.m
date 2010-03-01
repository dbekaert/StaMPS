function [ph_uw,msd]=uw_3d(ph,xy,day,ifgday_ix,options)
%UW_3D unwrap phase time series 
%
%   Andy Hooper, Jun 2007

% ============================================================
% 08/2009 AH: Goldstein alpha value added to options
% ============================================================

if nargin<3
    help uw_3d
    error('not enough arguments')
end

if nargin<4
    ifgday_ix=[];
end

if ~isfield(options,'master_day')
    options.master_day=0;
end

if ~isfield(options,'grid_size')
    options.grid_size=200;
end

if ~isfield(options,'prefilt_win')
    options.prefilt_win=32;
end

if ~isfield(options,'time_win')
    options.time_win=180;
end

if ~isfield(options,'unwrap_method')
    options.unwrap_method='3D';
end

if ~isfield(options,'goldfilt_flag')
    options.goldfilt_flag='y';
end

if ~isfield(options,'lowfilt_flag')
    options.lowfilt_flag='y';
end

if ~isfield(options,'gold_alpha')
    options.gold_alpha=0.8;
end

uw_grid_wrapped(ph,xy,options.grid_size,options.prefilt_win,options.goldfilt_flag,options.lowfilt_flag,options.gold_alpha);
uw_interp;
if isempty(ifgday_ix)
    uw_unwrap_space_time(day,options.unwrap_method,options.time_win,options.master_day);
else
    uw_sb_unwrap_space_time(day,ifgday_ix,options.unwrap_method,options.time_win);
end
uw_stat_costs(options.unwrap_method);
[ph_uw,msd]=uw_unwrap_from_grid(ph,xy,options.grid_size);

