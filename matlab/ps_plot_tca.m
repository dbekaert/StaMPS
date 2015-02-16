function [aps_corr,fig_name_tca,aps_flag] = ps_plot_tca(aps,aps_flag)
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
% 08/2014   DB      Include MODIS recalibrated support


if ischar(aps_flag)
   if strcmp(aps_flag,'a_l')==1
       % aps topo correlated linear correction
       aps_flag=1;
   elseif strcmp(aps_flag,'a_p')==1
       % aps topo correlated powerlaw correction
       aps_flag=2;
   elseif strcmp(aps_flag,'a_m')==1
       % aps topo correlated meris correction
       aps_flag=3;
   elseif strcmp(aps_flag,'a_e')==1
       % aps topo correlated ERA-I  correction
       aps_flag=4;
   elseif strcmp(aps_flag,'a_eh')==1
       % aps hydrostatic ERA-I correction
       aps_flag=5;
   elseif strcmp(aps_flag,'a_ew')==1
       % aps topo correlated ERA-I correction
       aps_flag=6;
   elseif strcmp(aps_flag,'a_w')==1
       % aps topo correlated WRF correction
       aps_flag=7;
   elseif strcmp(aps_flag,'a_wh')==1
       % aps hydrostatic WRF correction
       aps_flag=8;
   elseif strcmp(aps_flag,'a_ww')==1
       % aps topo correlated WRF correction
       aps_flag=9;
   elseif strcmp(aps_flag,'a_mi')==1
       % aps topo correlated  MERIS (non-interpolated)
       aps_flag=10;
   elseif strcmp(aps_flag,'a_pk')==1
       % Spatial maps of the coefficient relating phase and tropo for power-law
       aps_flag=11;
   elseif strcmp(aps_flag,'a_M')==1
       % aps topo correlated modis correction
       aps_flag=12;
   elseif strcmp(aps_flag,'a_MI')==1
       % aps topo correlated modis (non-interpolated)
       aps_flag=13;       
   elseif strcmp(aps_flag,'a_m+a_eh')==1
       % aps topo correlated MERIS plus a hydrostatic component from ERA-I
       aps_flag=14;
   elseif strcmp(aps_flag,'a_mi+a_eh')==1
       % aps topo correlated MERIS (non-interpolated) plus a hydrostatic component from ERA-I
       aps_flag=15;
   elseif strcmp(aps_flag,'a_M+a_eh')==1
       % aps topo correlated modis plus a hydrostatic component from ERA-I
       aps_flag=16;
   elseif strcmp(aps_flag,'a_MI+a_eh')==1
       % aps topo correlated modis (non-interpolated) plus a hydrostatic component from ERA-I
       aps_flag=17;
   elseif strcmp(aps_flag,'a_lman')==1
        % aps topo correlated manually estimated
       aps_flag=18;
   elseif strcmp(aps_flag,'a_RM')==1
       % aps topo correlated modis recalibrated correction
       aps_flag=19;
   elseif strcmp(aps_flag,'a_RMI')==1
       % aps topo correlated modis recalibrated (non-interpolated)
       aps_flag=20;          
   elseif strcmp(aps_flag,'a_RM+a_eh')==1
       % aps topo correlated modis recalibrated plus a hydrostatic component from ERA-I
       aps_flag=21;
   elseif strcmp(aps_flag,'a_RMI+a_eh')==1
       % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from ERA-I
       aps_flag=22;       
   end
end




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
elseif aps_flag==18 % current implementation of aps correction (manually estimated)
    aps_corr = aps.strat_corr;
    fig_name_tca = ' (linear)';    
elseif aps_flag==19 % modis correction
    aps_corr = aps.ph_tropo_modis_recal;
    fig_name_tca = ' (modis recal)';
elseif aps_flag==20 % modis correction (not interpolated)
    aps_corr = aps.ph_tropo_modis_no_interp_recal;
    fig_name_tca = ' (modis recal)';   
elseif aps_flag==21 % modis correction + ERA hydro component
    ix_no_modis = find(sum(aps.ph_tropo_modis)==0);
    aps_corr = aps.ph_tropo_modis_recal+aps.ph_tropo_era_hydro;
    aps_corr(:,ix_no_modis)=0;
    fig_name_tca = ' (modis recal + ERA hydro)';    
elseif aps_flag==22 % modis correction (not interpolated)
    aps_corr = aps.ph_tropo_modis_no_interp_recal+aps.ph_tropo_era_hydro;
    fig_name_tca = ' (modis recal + ERA hydro)';     
else
    error('not a valid APS option')
end


