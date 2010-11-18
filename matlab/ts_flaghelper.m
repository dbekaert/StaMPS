% TS_flaghelper for PS_PLOT function
%
%
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   11/2010 MMC & MA: plotting time series using 'ts' option 
%   ======================================================================


savename=['ps_plot_ts_',value_type];
savetxt='ps_plot_ts_matname.txt';  % if you change here update ts_plot as well
% save required matrices for TS PLOT
% save ps_ts ph_uw unwrap_ifg_index ref_ps day n_ps lambda ifg_list master_day m ph_all G
% put matrice name to 
dlmwrite(savetxt,[savename, '.mat'],'') % put savename to ts_plot_matname.txt
save(savename,'ph_uw','lonlat', 'unwrap_ifg_index', 'ref_ps',...
    'day', 'n_ps', 'lambda', 'ifg_list', 'master_day') % 'ph_all', 'G', 'm'
%clear all % clean up 

%EOF