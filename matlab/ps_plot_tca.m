function [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag)
% [aps_corr,fig_name_tca] = ps_plot_tca(aps,aps_flag)
% function that figes the selected aps correction and a string describing 
% the correction
%
% By David Bekaert - University of leeds
% October 2013
% modifications:
% 11/2013   DB      Include WRF and update naming dry to hydro
% 11/2013   DB      Include the non-interpolated version of MERIS
% 02/2014   DB      Include the plotting of the powerlaw spatial maps of K
% 05/2014   DB 	    Include MODIS support


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
    aps_corr = aps.ph_tropo_era_hydro;
    fig_name_tca = ' (era hydro)';
elseif aps_flag==6 % ERA-I correction
    aps_corr = aps.ph_tropo_era_wet;
    fig_name_tca = ' (era wet)';
elseif aps_flag==7 % WRF correction
    aps_corr = aps.ph_tropo_wrf;
    fig_name_tca = ' (wrf)';
elseif aps_flag==8 % ERA-I correction
    aps_corr = aps.ph_tropo_wrf_hydro;
    fig_name_tca = ' (wrf hydro)';
elseif aps_flag==9 % ERA-I correction
    aps_corr = aps.ph_tropo_wrf_wet;
    fig_name_tca = ' (wrf wet)';
elseif aps_flag==10 % meris correction (not interpolated)
    aps_corr = aps.ph_tropo_meris_no_interp;
    fig_name_tca = ' (meris)';
elseif aps_flag==11 % power-law spatial maps of K
    aps_corr = aps.K_tropo_powerlaw;
    fig_name_tca = ' (powerlaw - spatial K map)';
elseif aps_flag==12 % modis correction
    aps_corr = aps.ph_tropo_modis;
    fig_name_tca = ' (modis)';
elseif aps_flag==13 % modis correction (not interpolated)
    aps_corr = aps.ph_tropo_modis_no_interp;
    fig_name_tca = ' (modis)';
elseif aps_flag==14 % meris correction + ERA hydro component
    ix_no_meris = find(sum(aps.ph_tropo_meris)==0);
    aps_corr = aps.ph_tropo_meris+aps.ph_tropo_era_hydro;
    aps_corr(:,ix_no_meris)=0;
    fig_name_tca = ' (meris + ERA hydro)';
elseif aps_flag==15 % meris correction (not interpolated) + ERA hydro component
    aps_corr = aps.ph_tropo_meris_no_interp + aps.ph_tropo_era_hydro;
    fig_name_tca = ' (meris + ERA hydro)';
elseif aps_flag==16 % modis correction + ERA hydro component
    ix_no_modis = find(sum(aps.ph_tropo_modis)==0);
    aps_corr = aps.ph_tropo_modis+aps.ph_tropo_era_hydro;
    aps_corr(:,ix_no_modis)=0;
    fig_name_tca = ' (modis + ERA hydro)';    
elseif aps_flag==17 % modis correction (not interpolated)
    aps_corr = aps.ph_tropo_modis_no_interp+aps.ph_tropo_era_hydro;
    fig_name_tca = ' (modis + ERA hydro)';    
else% current implementation of aps correction
    aps_corr = aps.strat_corr;
    fig_name_tca = ' (linear)';
end


