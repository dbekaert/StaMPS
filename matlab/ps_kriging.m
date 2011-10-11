function [sol var_sol min_dist mean_dist nof_obs] = ps_kriging(x,x0,Nmax,dx_max,kriging_method,cv_model)

% Kriging algorithm (for general use)
%
% Input:    - x               data points [cn1 cn2 z]
%           - x0              new points [cn1 cn2]
%           - Nmax            number of nearest neighbors to use
%           - dx_max          search radius for neighbors
%           - kriging_method  kriging method,
%                             1 = Simple Kriging
%                             2 = Ordinary Kriging with one nonbias
%                             constraint
%                             3 = Ordinary Kriging with number of
%                             parameters nonbias constraints
%                             4 = Universal Kriging with drift of
%                             order 1
%                             5 = Universal Kriging with drift of
%                             order 2
%           - cv_model        covariance model
%                             1 = nugget effect
%                             2 = exponential model
%                             3 = gaussian model
%                             4 = spherical model
%                             (new models can easily be added)
%									   e.g.  cv_model = [1 NaN c2;covariance_model a c1];
% Output:   - sol             solution
%
% ----------------------------------------------------------------------
% File............: ps_kriging.m
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


% -----------------------------------------------------
% Initialize
% -----------------------------------------------------
%keyboard
Nmodel = size(cv_model,1);
Npar = 1;
Nx0 = size(x0,1);
Nx = size(x,1);
Nmax = min(Nmax,Nx);

% -----------------------------------------------------
% Select data points
% -----------------------------------------------------
if dx_max~=0
  if size(x0,1)>1
	 x0_mean = mean(x0,1);
  else
   x0_mean=x0;
  end
	dx = sqrt((x(:,1)-x0_mean(1)).^2+(x(:,2)-x0_mean(2)).^2);
	[dx_sort,index] = sort(dx);
  %keyboard
	index2 = find(dx_sort(1:Nmax)<dx_max);
	index = index(index2);
	x = x(index,:);
	Nx = size(x,1);
end


if Nx<3
  kriging_method=2; %change to ordinary kriging, otherwise singular
                    %k matrix
end

% -----------------------------------------------------
% Setup Kriging system
% -----------------------------------------------------
xa = repmat(x(:,1),1,Nx);
xb = xa';
ya = repmat(x(:,2),1,Nx);
yb = ya';
x0a = repmat(x(:,1),1,Nx0);
x0b = repmat(x0(:,1)',Nx,1);
y0a = repmat(x(:,2),1,Nx0);
y0b = repmat(x0(:,2)',Nx,1);
%distances between observations
h = double(sqrt((xa-xb).^2+(ya-yb).^2));
%distances from obs to point of interest
h0 = double(sqrt((x0a-x0b).^2+(y0a-y0b).^2));

min_dist=min(h0);
mean_dist=mean(h0);
nof_obs=length(find(h0<cv_model(2,2)));
%This is the redundacy (Cm)  matrix contaninig the covariances of the observations
%Only the diagonal includes the nugget because the non-diagonal are =zero for the nugget variance function
k = zeros(Nx,Nx);
%This is the proximity vector (d), contaninig the initial weight
%assigned to each observation obtained from the covariance function evaluated at the wanted location.
%The nugget is not included in k0 to filter out the noise component.
k0 = zeros(Nx,Nx0);

for w = 1:Nmodel
  k = k+covar_function(h,cv_model(w,:));
%  keyboard
  if cv_model(w,1)~=1
    k0 = k0+covar_function(h0,cv_model(w,:));
  end	
end


switch kriging_method
  case 1 %Simple Kriging
    avg = mean(x(:,3));
    x(:,3) = x(:,3)-avg;

  case 2 %Ordinary Kriging (one nonbias constraint)
    t = ones(Nx*Npar,1);
    k = [k t;t' 0];
    t = ones(1,Nx0*Npar);
    k0 = [k0;t];

    avg = mean(x(:,3));
    x(:,3) = x(:,3)-avg;
    
  case 3 %Ordinary Kriging (Npar nonbias constraints)
    error(['Sorry, Kriging with multiple parameters is not ' ...
           'implemented yet']);
    t = kron(ones(Nx,1),eye(Npar));
    k = [k t;t' zeros(Npar)];
    t = kron(ones(1,Nx0),eye(Npar));
    k0 = [k0;t];
    
  case 4 %Universal Kriging with drift of order 1
    t = kron(ones(Nx,1),eye(Npar));
    t2 = kron(x(:,1:2),eye(Npar));
    k = [k t t2; [t' zeros(Npar,3*Npar)]; [t2' zeros(2*Npar,3*Npar)]];
    t = kron(ones(1,Nx0),eye(Npar));
    t2 = kron(x0',eye(Npar));
    k0 = [k0;t;t2];
    
  case 5 %Universal Kriging with drift of order 2
    error(['Sorry, Universal Kriging with drift of order 2 is ' ...
           'not implemented yet']);
    
  otherwise
    error('You specified a wrong Kriging method');
end


% -----------------------------------------------------
% Solve the Kriging system
% -----------------------------------------------------

%l = inv(k)*k0;
%keyboard
l=k\k0; %a little bit faster (mcc)
%l are the weights of the observations when evaluated at the needed location (p0)
%sol(p0)= sum( weights(p0)*data);

%mcc
sol = l(1:Nx,:)'*x(:,3);
var_sol=repmat(0,size(sol));
%The error variance var_sol= C00 - sum(wi*Ci0) - lambda (lagrange param)
%the same as var_sol= C00 - sum (wi*di)
%because the last term of di=1 and wi=lambda

for n=1:length(sol)
	var_sol(n)=k(1,1)-l(:,n)'*k0(:,n);%where k= Cn and k(1,1)=C00; l represent the weihgts and k0 the proximity vector
end

switch kriging_method
  case 1 %Simple Kriging
    sol = sol+avg;
  case 2 %Ordinary Kriging (one nonbias constraint)
    sol = sol+avg;
end

if isnan(sol)
 var_sol=NaN;
 min_dist=NaN;
 mean_dist=NaN;
 warning('Solution is NaN');
 
end
%subfunction
function cv = covar_function(h,cv_model)

model = cv_model(1);
a = cv_model(2);
c = cv_model(3);

switch model
  case 1 %nugget
    cv = kron(eye(size(h)),c);
  case 2 %exponential
    cv = c*exp(-h/a);
  case 3 %gaussian
    cv = c*exp(-(h/a).^2);
  case 4 %spherical
    cv = c*(1-1.5*min(h/a,1)+0.5*min(h/a,1).^3);
  otherwise
    error('You specified a wrong covariance model');
end

