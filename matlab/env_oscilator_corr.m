function [oscilatior_corr_ifgs,oscilatior_corr_velocity] = env_oscilator_corr(envisat_flag,forced_sm_flag)
% This function perform oscialtor drift correction for envisat interferograms based on
% Petar Marinkovic presentation at ESA Living Planet 2013 in Edinburgh.
% Approximation formula:
% dR/year = c/2 * (slantRangeTime_FAR - slantRangeTime_NEAR) * corrPerYear
% 'Apparent displacement' correction:
% dR/year ~ (7.8m * 5000)*3.87e-7 ~ 0.01482m
% Correction uses the pixel information in range and in adition for
% interferograms the temporal baseline information
%
% OPTIONAL INPUT:
% envisat_flag                  'y' when envisat or 'n'. When empty, a
%                               search is performed in the master.res file
%
% OUTPUTS:
% oscilatior_corr_velocity      Correction in mm for the velocity
% oscilatior_corr_ifgs          Correction for individual interferograms in rad
%
% Correction is defined such that:
% Corrected interferogram/velocity = original interferogram/velocity - correction
%
% P. Marinkovic Envisat oscialtor drift correction coded by David Bekaert -- University of Leeds 2014
% 
% cite as:
% P. Marinkovic (PPO.labs) and Y. Larsen (NORUT)
% Consequences of Long-Term ASAR Local Oscillator Frequency Decay - an Empirical Study of 10 Years of Data 
% ESA Living Planet Symposium (2013)
%
% Modifications
% 04/2014       DB      Put non-envisat fix to zeros
% 04/2014       DB      Allow forced SM oscialtor drift computation
% 06/2014       DB      Fixed error for SM computation and added extra envisat check.
% 06/2014       DB      Fix in case not envisat and forced SM.
% 11/2014       DB      Fix to make windows and linux compatible 
% 03/2015       DB      Clean script output
% 01/2017       DB      Bug fix for non-envisat SM case. n_ifg was not n_image


if nargin<1 || isempty(envisat_flag)
  % checking if this is envisat or not
  platform=getparm('platform');
  if isempty(platform)
    if exist('master.res','file')==2
        master_file = 'master.res';
    elseif exist('../master.res','file')==2
        master_file = '../master.res';
    else
        master_file = [];
    end

    if ~isempty(master_file)
        
        % make windows and linux compatible
        a = fileread(master_file);
        ix = strfind(a,'ASAR');
       if ~isempty(ix)
           platform='ENVISAT';
       end
    else
        fprintf('Could not check if this is Envisat \n')
    end
  end
  if strcmpi(platform,'ENVISAT')
      envisat_flag = 'y';
      fprintf('This is Envisat, oscilator drift is being removed... \n')
  else
       envisat_flag = 'n';
  end

end

small_baseline_flag = getparm('small_baseline_flag');
if nargin<2 
    forced_sm_flag=0;
end

if forced_sm_flag==1
    small_baseline_flag='n';
end



if strcmp(envisat_flag,'y')
    load psver
    % use the ps2.mat data
    ps = load(['ps' num2str(psver) '.mat']);
    lambda = getparm('lambda');

    
    % velocity map correction:
    envisat_resolution = 7.8;                   % ground range resolution [m]
    Oscilator_drift_corr_year = 3.87e-7;         % drift correction in range [1/year]      

    % velocity correction in mm
    oscilatior_corr_velocity = (envisat_resolution*ps.ij(:,3))*Oscilator_drift_corr_year*1000;

    % interferogram
    if strcmp(small_baseline_flag,'y')
        n_ifg = ps.n_ifg;
        delta_year = (ps.ifgday(:,2)-ps.ifgday(:,1))./365.25;
    else
        n_ifg = ps.n_image;
        delta_year = (ps.day-ps.master_day)./365.25;
    end
    oscilatior_corr_ifgs = -4.*pi./lambda.*repmat(oscilatior_corr_velocity,1,n_ifg)./1000.*repmat(delta_year',ps.n_ps,1);

    
else
    load psver
    ps = load(['ps' num2str(psver) '.mat']);
    
    if strcmp(small_baseline_flag,'y')
        oscilatior_corr_ifgs = zeros([ps.n_ps ps.n_ifg]);
    else
        n_ifg = ps.n_image;     % bug fix DB
        oscilatior_corr_ifgs = zeros([ps.n_ps n_ifg]);
    end
    oscilatior_corr_velocity = zeros([ps.n_ps 1]);

end