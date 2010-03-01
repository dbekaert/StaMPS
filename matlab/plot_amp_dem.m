function []=plot_amp_dem(az_down,rg_right,exponent,scale_factor,use_slope)
%PLOT_AMP_DEM combine amp and dem and display
%    MAKE_AMP_DEM(BLUE_DOWN,BLUE_RIGHT,RED_CONTRAST,BLUE_BRIGHTNESS)
%    RED_CONTRAST  (default 0.5) = contrast of interferogram amplitude
%    BLUE_BRIGHTNESS (default 1) = brightness of DEM or DEM slope
%    USE_SLOPE (default 1) = set to 0 to plot DEM elevation instead 
%
%   Andy Hooper, Dec 2007
%
% ======================================================
%   10/2008 AH: Generalise for other sensors
% ======================================================


if nargin<1
    az_down=0;
end

if nargin<2
    rg_right=0;
end

if nargin<3
    exponent=0.5;
end

if nargin<4
    scale_factor=1;
end

if nargin<5
    use_slope=1;
end

[dd,ar]=combine_amp_dem(az_down,rg_right,exponent,scale_factor,use_slope);

image(dd)

load sar_parms.in -ascii
prf=sar_parms(1);
rsr=sar_parms(2)*1e6;

%disp(['If looks OK, add ',num2str(-(az_down)*5/prf),' to azimuth time in master.res'])
if rg_right ~=0 || az_down ~=0
   disp(['If looks OK, increment values in dem.dorisin and geocode.dorisin by: '])
   disp(['M_RG_T_ERROR    ',num2str(-(rg_right-0.5)/rsr/2,6)])
   disp(['M_AZ_T_ERROR    ',num2str(-(az_down)*ar/prf,6)])
end
