function [q]=K2q(K,range_pixel,sensor)
%Convert from K (slope of phase with perpendicular baseline) to DEM error
%   q=K2q(K,range_pixel,sensor)
%   where K is in rad/m, range_pixel refers to pixel number in range with no 
%   oversampling and sensor is 'E' for ERS/Envisat I2, 'A' for ALOS FBS 
%   or 'R2U18' for Radarsat-2 ultrafine beam 18. 
%   q the DEM error in m.
%
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   11/2010 DB: ALOS option added
%   11/2011 AH: Radarsat-2 Ultrafine 18 option added
%   01/2011 AH: Sign changed so that DEM error is +ve if reality > DEM
%   ======================================================================

if nargin<2
    range_pixel=2300;
end
if nargin<3
    sensor='E';
end

if strcmpi(sensor,'A')
    %ALOS
    slant_range=4.684;
    near_range=846014.31;
    h=691650;
    lambda= 0.236057;
elseif strcmpi(sensor,'R2U18')
    %RADARSAT-2 ultrafine 18
    slant_range=1.330719;
    near_range=1042559.7347;
    h=786503;
    lambda =0.05546577;
else
    %ERS/ENVISAT
    slant_range=7.9048902811596;
    near_range=831265;
    h=794141;
    lambda =0.0565647;
end

rho = near_range+range_pixel*slant_range;

re=6345245;
%load parms

alpha=pi-acos((rho.^2+re^2-(re+h)^2)/2./rho/re);

q=-K*lambda.*rho.*sin(alpha)/4/pi;
