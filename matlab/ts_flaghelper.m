% TS_flaghelper for PS_PLOT function
%
%
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   11/2010 MMC & MA: plotting time series using 'ts' option 
%   02/2015 AH Store phase in mm instead of radians
%   ======================================================================


savename=['ps_plot_ts_',value_type];
savetxt='ps_plot_ts_matname.txt';  % if you change here update ts_plot as well
% save required matrices for TS PLOT
% save ps_ts ph_uw unwrap_ifg_index ref_ps day n_ps lambda ifg_list master_day m ph_all G
% put matrice name to 
dlmwrite(savetxt,[savename, '.mat'],'') % put savename to ts_plot_matname.txt
load ps2 bperp;
lambda=getparm('lambda');
ph_mm=-ph_uw*lambda*1000/(4*pi);
save(savename,'ph_mm','lonlat', 'unwrap_ifg_index', 'ref_ps',...
    'day', 'n_ps', 'lambda', 'ifg_list', 'master_day', 'bperp') % 'ph_all', 'G', 'm'
%clear all % clean up 

%EOF