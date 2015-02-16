function [aps_corr,fig_name_tca,iono_flag] = ps_plot_ica(aps,iono_flag)
% [aps_corr,fig_name_tca] = ps_plot_ica(aps,iono_flag)
% function that gives the selected ionopsheric correction and a string describing 
% the correction
%
% By David Bekaert - University of leeds
% October 2014
% modifications:



if ischar(iono_flag)
   if strcmp(iono_flag,'i_az')==1
       % azimuth shift method for ionospheric delays
       iono_flag=1;
   end
end



if iono_flag==1 % linear correction
    aps_corr = aps.ph_iono_azshift;
    fig_name_tca = ' (azimuth shift method)';
else
    error('not a valid Ionospheric option')
end


