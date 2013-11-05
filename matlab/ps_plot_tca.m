function [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag)
% [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag)
% function that figes the selected aps correction and a string describing 
% the correction
%
% By David Bekaert - University of leeds
% October 2013

if aps_flag==1 % linear correction
    aps_corr = aps.ph_tropo_linear;
    fig_name_tca = ' (linear)';
elseif aps_flag==2 % powerlaw correlation
    aps_corr = aps.ph_tropo_powerlaw;
    fig_name_tca = ' (powerlaw)';
elseif aps_flag==3 % meris correction
    aps_corr = aps.ph_tropo_meris;
    fig_name_tca = ' (meris)';
elseif aps_flag==4 % ERA-I correction
    aps_corr = aps.ph_tropo_era;
    fig_name_tca = ' (era)';
elseif aps_flag==5 % ERA-I correction
    aps_corr = aps.ph_tropo_era_dry;
    fig_name_tca = ' (era dry)';
elseif aps_flag==6 % ERA-I correction
    aps_corr = aps.ph_tropo_era_wet;
    fig_name_tca = ' (era wet)';
else% current implementation of aps correction
    aps_corr = aps.strat_corr;
    fig_name_tca = '';
end


