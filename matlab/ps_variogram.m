function [y,dyda] = ps_variogram(model,h,a,b,c1,c2)

% Function for the evaluation of variogram model (used by
% ps_mrqmin.m, Levenberg-Marquardt non-linear LS)
%
% Input:    
%
% Output:   
%
% ----------------------------------------------------------------------
% File............: ps_variogram.m
% Version & Date..: 1.7.1, 01-NOV-2006
% Authors.........: Freek van Leijen
%                   Delft Institute of Earth Observation and Space Systems
%                   Delft University of Technology
% ----------------------------------------------------------------------
%
% This software is developed by Delft University of Technology and is
% intended for scientific use only. Applications for commercial use are
% prohibited.
%
% Copyright (c) 2004-2006 Delft University of Technology, The Netherlands
%

Nh = length(h);
if Nh>1&model==4
  error(['For now only implemented for 1 variable (due to spherical ' ...
         'model).']);
end
dyda = repmat(NaN,Nh,4);

if model==5
model=2
end
switch model
  case 1
    error(['You can not estimate the nugget independently, choose ' ...
           'another model (nugget included)']);
    
  case 2 %exponential
    y = c1*(1-exp(-h/a))+c2;
    dyda(:,1) = -c1*h/a^2 .*exp(-h/a);
    dyda(:,2) = 0;   
    dyda(:,3) = y/c1;
    dyda(:,4) = 0;
    
  case 3 %gaussian
    y = c1*(1-exp(-(h/a).^2))+c2;
    dyda(:,1) = -2*c1*h.^2/a^3 .* exp(-(h/a).^2);
    dyda(:,2) = 0;
    dyda(:,3) = y/c1;
    dyda(:,4) = 0;
  
  case 4 %spherical
    y = c1*(1.5*min(h/a,1)-0.5*min(h/a,1).^3)+c2;
    if h<a
      dyda(:,1) = c1*(-1.5*h/a^2+1.5*h.^3/a^4);
      dyda(:,2) = 0;
      dyda(:,3) = y/c1;
      dyda(:,4) = 0;
    else
      dyda(:,1) = 0;
      dyda(:,2) = 0;
      dyda(:,3) = 0;
      dyda(:,4) = 0;
    end
  
  otherwise
    error('You specified a wrong variogram model');
    
end

