function []=plot_sb_baselines(vargin)
% PLOT_SB_BASELINES replaced by SB_BASELINE_PLOT
% Optional an input argument ix can be specified containing a vector of the
% small baseline interferograms to keep which will be plotted in the baseline plot.
%
% Modifications:
% DB 12/2012    fix error for input arguments and add the some syntax
% DB 01/2016    Suppress the warning. Alternative remove the function.

% warning('name of this function changed to sb_baseline_plot')

if nargin>0                     % [DB] changed to nargin > 0, as ix cannot be specified otherwize
    sb_baseline_plot(vargin)
else
    sb_baseline_plot
end
