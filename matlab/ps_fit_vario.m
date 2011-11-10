function [a,b,c1,c2,emp_vario] = ps_fit_vario(x,y,z,model,a0,b0,c10,c20,decor_dist,detail_plots,Nlags,previous_vario,name_plot)

% Estimation of variogram model
%
% Input:    - x                 x-coordinates
%           - y                 y-coordinates
%           - z                 z-values
%           - model             variogram model:
%                               2 = exponential
%                               3 = gaussian
%                               4 = spherical
%           - a0                initial value range
%           - b0                initial value hole effect
%           - c10               initial value sill
%           - c20               initial value nugget
%
% Output:   - a                 estimated range
%           - b                 estimated hole effect
%           - c1                estimated sill
%           - c2                estimated nugget
%           - emp_vario         empirical variogram
%
% ----------------------------------------------------------------------
% File............: ps_fit_vario.m
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



% -----------------------------------------------------
% Initialize
% -----------------------------------------------------

%global fig detail_plots visible_plots
if nargin<10

	detail_plots='n';

end
c10=double(c10);
c20=double(c20);
% -----------------------------------------------------
% Calculate raw variogram
% -----------------------------------------------------

Nx = size(x,1);
%xa = repmat(x,1,Nx);
%xb = xa';
%ya = repmat(y,1,Nx);
%yb = ya';
%za = repmat(z,1,Nx);
%zb = za';
%memory efficient
%h1=speye(Nx);
%h2=speye(Nx);

h1=repmat(single(0),Nx,Nx);
h2=repmat(single(0),Nx,Nx);

for n =1:Nx
 
	h1(n,:)=sqrt((x-x(n)).^2+(y-y(n)).^2);
  %ind_out=find(h1(n,:)>20000);
 % h1(n,ind_out)=zeros(size(ind_out));
  
	h2(n,:)=0.5*((z-z(n)).^2);
 % h2(n,ind_out)=zeros(size(ind_out));
end

%h1 = sqrt((xa-xb).^2+(ya-yb).^2);
%h2 = 0.5*((za-zb).^2);

%h1 = h1(:);
%h2 = h2(:);

h1(1:Nx+1:Nx*Nx) = []; %remove h1==0
h2(1:Nx+1:Nx*Nx) = [];

%if strcmp(detail_plots,'y')
%  fig = fig+1;
%  figure(fig);hold on
%  if strcmp(visible_plots,'n')
%    set(gcf,'visible','off');
%  end
%  plot(h1,h2,'*');
%  xlabel('Distance [m]')
%  ylabel('Variogram [rad^2]')
%end


% -----------------------------------------------------
% Calculate experimental variogram
% -----------------------------------------------------
if ~exist('decor_dist')
	max_dist = max(h1(:));
else
	max_dist=decor_dist;

end

if ~exist('Nlags')
lags = [0:max_dist/50:max_dist];
Nlags = length(lags)-1;

else
lags = [0:max_dist/(Nlags-1):max_dist];
Nlags=Nlags-1;
end

emp_vario = repmat(NaN,Nlags,3);
for v = 1:Nlags
    index = find((h1>=lags(v))&(h1<lags(v+1)));
    emp_vario(v,1) = length(index);
    if emp_vario(v,1)~=0
      emp_vario(v,2) = nanmean(h2(index));
    end
    emp_vario(v,3) = (lags(v)+lags(v+1))/2;
end
index = find(isnan(emp_vario(:,2)));
emp_vario(index,:) = [];
clear h1 h2
%keyboard
%if strcmp(detail_plots,'y')
%  figure;
%  plot(emp_vario(:,3),emp_vario(:,2),'r*','markersize',10)
%	print('-dpng','emp_vario.png');
%end

if size(emp_vario,1)<round(Nlags/1.75)
	last_lag =  size(emp_vario,1);

else
	last_lag = round(Nlags/1.75);

end

emp_vario = emp_vario(1:last_lag,:);
Nlags = size(emp_vario,2);

% -----------------------------------------------------
% Estimate experimental variogram
% -----------------------------------------------------

lam = [a0 c10 c20];

if exist('decor_dist')
	ind_max=max(find(emp_vario(:,3)<decor_dist));
else
	ind_max=size(emp_vario,1);
end
if exist('previous_vario')
  if ~isempty(previous_vario)
[y_final,dyda] = ps_variogram(previous_vario(1),emp_vario(:,3),previous_vario(2),previous_vario(3),previous_vario(4),previous_vario(5));
emp_vario(:,2)=emp_vario(:,2)-y_final;
  end
end

switch model
  case 1
    error(['You can not estimate the nugget independently, choose ' ...
           'another model (nugget included)']);
  case 2 %exponential                %lsqcurvefit uses a function handle (like a pointer to a function given by @) to estimate th best function parameters that fit the data for a starting set of parameters (given in lam)
    [lam,resnorm,residual,exitflag] = lsqcurvefit(@ps_variogram_exp,lam,emp_vario(1:ind_max,3),emp_vario(1:ind_max,2),[0 0 0],[],optimset('Display','off'));
  case 3 %gaussian
    [lam,resnorm,residual,exitflag] = lsqcurvefit(@ps_variogram_gaus,lam,emp_vario(1:ind_max,3),emp_vario(1:ind_max,2),[0 0 0],[],optimset('Display','off'));
  case 4 %spherical
    [lam,resnorm,residual,exitflag] = lsqcurvefit(@ps_variogram_spher,lam,emp_vario(1:ind_max,3),emp_vario(1:ind_max,2),[0 0 0],[],optimset('Display','off'));

  otherwise
    error('You specified a wrong variogram model');
end
    
%if exitflag==0
%  warning('Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.');
%elseif exitflag==-1
%  warning('Algorithm was terminated by the output function.');
%elseif exitflag==-2
%  warning('Problem is infeasible: the bounds lb and ub are inconsistent.');
%elseif exitflag==-3
%  warning('Optimization could not make further progress.');
%end

a = lam(1);
b = 0;
c1 = lam(2);
c2 = lam(3);


if strcmp(detail_plots,'y')
  [y_final,dyda] = ps_variogram(model,emp_vario(:,3),a,b,c1,c2);
  figure;
  plot(emp_vario(:,3),y_final,'r','linewidth',2);
  hold on
  plot(emp_vario(:,3),emp_vario(:,2),'r*')
	hold off
  if exist('name_plot')
   print('-dpng',['emp_vario_'  name_plot  '.png']);
  else
   print('-dpng','emp_vario.png');
  end
end
%keyboard


% subfunctions



function y = ps_variogram_exp(x,xdata);

y = x(2)*(1-exp(-xdata/x(1)))+x(3);


function y = ps_variogram_gaus(x,xdata);

y = x(2)*(1-exp(-(xdata/x(1)).^2))+x(3);


function y = ps_variogram_spher(x,xdata);

y = x(2)*(1.5*min(xdata/x(1),1)-0.5*min(xdata/x(1),1).^3)+x(3);
