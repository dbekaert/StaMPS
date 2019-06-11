%   parse parameters for PS_PLOT function
%
%   Andy Hooper, June 2006
%
%   modifications:
%   11/2010 MA: Initial support for 'ts' option.
%   04/2013 DB: Include support of topocorrelated aps 'a_m' for merris, 
%               'a_l' for linear, and 'a_p' for powerlaw option.
%   04/2013 DB: Include bandfilter option 'ifg i'.
%   05/2013 DB: Allow ext data path to be a folder or a file.
%   06/2013 DB: Allows units to be specified 
%   02/2014 DB: Allow for spatial maps of K to be plotted for the power-law technique
%   05/2014 DB: Include MODIS support
%   08/2014 DB: MODIS recalibrated support
%   10/2014 DB: Support for ionospheric delays
%   03/2014 DB: Fix for a bug introduced when APS option does not exist.
%               Add the missing option too  
%   08/2016 AH: Remove "break" that is incompatible with Matlab2015
%   11/2017 DB: Change the APS naming convention to allow more models, add
%               merra and merra2 models
%   11/2017 DB: Adding the GACOS model 
%   02/2018 DB: Adding the NARR model and adding in backward compatibility


% list of arguments for ps_plot excluding value type
arglist={'plot_flag','lims','ref_ifg','ifg_list','n_x','cbar_flag',...
    'textsize','textcolor','lon_rg','lat_rg','units'};  % excluding value_type

%varargin  % debug
%stdargin=nargin ; 
%fprintf('Number of inputs = %d\n', nargin)

Noptargin = length(varargin(:));  % treat extra options as varargin, for ps_plotTS('v-d',4,'ts') is 3

% Defaults, later below update with user input.
ts_flag=0;              % when 1 do a time-series display
aps_flag=0;             % when 0 default aps estimation 
                        % when 1 linear correction aps
                        % when 2 powerlaw correction aps
                        % when 3 meris correction aps
iono_flag = 0;      
aps_band_flag=0;        % when 1 plot the aps estimated band frequencies
ifg_number = 0;         % the ifg number for the aps_band_flag
ext_data_flag=0;        % when 1 there will be external data to be plotted.
ext_data_path = [];     % path to the external data to be displayed
plot_flag=1; 



% search for char parameter like 'ts','a_l','a_m','a_p', 'a_e', 'a_M', 'a_m+a_eh' , 'ifg i^th', 'ext PATH'
prmsrch=[];
for k=1:Noptargin,   
   if strcmp(varargin{k},'ts')==1
       % time series plot
       ts_flag=1;   
       prmsrch=logical([prmsrch strcmp(varargin{k},'ts') ]);
   elseif strcmp(varargin{k},'a_linear')==1 || strcmp(varargin{k},'a_l')==1
       % aps topo correlated linear correction
       aps_flag=1;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_linear') ]);
   elseif strcmp(varargin{k},'a_powerlaw')==1 || strcmp(varargin{k},'a_p')==1
       % aps topo correlated powerlaw correction
       aps_flag=2;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_powerlaw') ]);
   elseif strcmp(varargin{k},'a_meris')==1 ||  strcmp(varargin{k},'a_m')==1
       % aps topo correlated meris correction
       aps_flag=3;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_meris') ]);
   elseif strcmp(varargin{k},'a_erai')==1 || strcmp(varargin{k},'a_e')==1
       % aps topo correlated ERA-I  correction
       aps_flag=4;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_erai') ]);
   elseif strcmp(varargin{k},'a_erai-h')==1 || strcmp(varargin{k},'a_eh')==1
       % aps hydrostatic ERA-I correction
       aps_flag=5;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_erai-h') ]);
   elseif strcmp(varargin{k},'a_erai-w')==1 ||  strcmp(varargin{k},'a_ew')==1
       % aps topo correlated ERA-I correction
       aps_flag=6;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_erai-w') ]);
   elseif strcmp(varargin{k},'a_wrf')==1 ||  strcmp(varargin{k},'a_w')==1
       % aps topo correlated WRF correction
       aps_flag=7;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_wrf') ]);
   elseif strcmp(varargin{k},'a_wrf-h')==1 || strcmp(varargin{k},'a_wh')==1
       % aps hydrostatic WRF correction
       aps_flag=8;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_wrf-h') ]);
   elseif strcmp(varargin{k},'a_wrf-w')==1 ||  strcmp(varargin{k},'a_ww')==1
       % aps topo correlated WRF correction
       aps_flag=9;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_wrf-w') ]);
   elseif strcmp(varargin{k},'a_meris-ni')==1 ||  strcmp(varargin{k},'a_mi')==1
       % aps topo correlated  MERIS (non-interpolated)
       aps_flag=10;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_meris-ni') ]);     
  elseif strcmp(varargin{k},'a_powerlaw-k')==1 ||  strcmp(varargin{k},'a_pk')==1
       % Spatial maps of the coefficient relating phase and tropo for power-law
       aps_flag=11;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_powerlaw-k') ]);   
 elseif strcmp(varargin{k},'a_modis')==1 || strcmp(varargin{k},'a_M')==1
       % aps topo correlated modis correction
       aps_flag=12;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_modis') ]);       
 elseif strcmp(varargin{k},'a_modis-ni')==1 || strcmp(varargin{k},'a_MI')==1
       % aps topo correlated modis (non-interpolated)
       aps_flag=13;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_modis-ni') ]);
       
  elseif strcmp(varargin{k},'a_meris+a_erai-h')==1 ||  strcmp(varargin{k},'a_m+a_eh')==1
       % aps topo correlated MERIS plus a hydrostatic component from ERA-I
       aps_flag=14;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_meris+a_erai-h')]);     
  elseif strcmp(varargin{k},'a_meris-ni+a_erai-h')==1 || strcmp(varargin{k},'a_mi+a_eh')==1
       % aps topo correlated MERIS (non-interpolated) plus a hydrostatic component from ERA-I
       aps_flag=15;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_meris-ni+a_erai-h')]);     
  elseif strcmp(varargin{k},'a_modis+a_erai-h')==1 || strcmp(varargin{k},'a_M+a_eh')==1
       % aps topo correlated modis plus a hydrostatic component from ERA-I
       aps_flag=16;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_modis+a_erai-h')]);     
  elseif strcmp(varargin{k},'a_modis-ni+a_erai-h')==1 || strcmp(varargin{k},'a_MI+a_eh')==1
       % aps topo correlated modis (non-interpolated) plus a hydrostatic component from ERA-I
       aps_flag=17;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_modis-ni+a_erai-h')]);    
  elseif strcmp(varargin{k},'a_linear-man')==1 || strcmp(varargin{k},'a_lman')==1
       % aps topo correlated linear manual
       aps_flag=18;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_linear-man')]);    
 elseif strcmp(varargin{k},'a_recalmodis')==1 || strcmp(varargin{k},'a_RM')==1
       % aps topo correlated modis recalibrated correction
       aps_flag=19;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_recalmodis') ]);       
 elseif strcmp(varargin{k},'a_recalmodis-ni')==1 ||  strcmp(varargin{k},'a_RMI')==1
       % aps topo correlated modis recalibrated (non-interpolated)
       aps_flag=20;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_recalmodis-ni') ]);
  elseif strcmp(varargin{k},'a_recalmodis+a_erai-h')==1 || strcmp(varargin{k},'a_RM+a_eh')==1
       % aps topo correlated modis recalibrated plus a hydrostatic component from ERA-I
       aps_flag=21;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_recalmodis+a_erai-h')]);     
  elseif strcmp(varargin{k},'a_recalmodis-ni+a_erai-h')==1 || strcmp(varargin{k},'a_RMI+a_eh')==1
       % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from ERA-I
       aps_flag=22;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_recalmodis-ni+a_erai-h')]);    
  elseif strcmp(varargin{k},'a_meris+a_wrf-h')==1 || strcmp(varargin{k},'a_m+a_wh')==1
       % aps topo correlated MERIS plus a hydrostatic component from WRF
       aps_flag=23;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_meris+a_wrf-h')]);     
   elseif strcmp(varargin{k},'a_meris-ni+a_wrf-h')==1  || strcmp(varargin{k},'a_mi+a_wh')==1
       % aps topo correlated MERIS (non-interpolated) plus a hydrostatic component from WRF
       aps_flag=24;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_meris-ni+a_wrf-h')]);     
   elseif strcmp(varargin{k},'a_modis+a_wrf-h')==1 ||  strcmp(varargin{k},'a_M+a_wh')==1
       % aps topo correlated modis plus a hydrostatic component from  WRF
       aps_flag=25;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_modis+a_wrf-h')]);     
   elseif strcmp(varargin{k},'a_modis-ni+a_wrf-h')==1 || strcmp(varargin{k},'a_MI+a_wh')==1
       % aps topo correlated modis (non-interpolated) plus a hydrostatic component from WRF
       aps_flag=26;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_modis-ni+a_wrf-h')]);    
   elseif strcmp(varargin{k},'a_recalmodis+a_wrf-h')==1 || strcmp(varargin{k},'a_RM+a_wh')==1
       % aps topo correlated modis recalibrated plus a hydrostatic component from WRF
       aps_flag=27;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_recalmodis+a_wrf-h')]);     
   elseif strcmp(varargin{k},'a_recalmodis-ni+a_wrf-h')==1 || strcmp(varargin{k},'a_RMI+a_wh')==1
       % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from WRF
       aps_flag=28;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_recalmodis-ni+a_wrf-h')]);  
  elseif strcmp(varargin{k},'a_merra')==1
        % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from WRF
        aps_flag=29;
        prmsrch=logical([prmsrch strcmp(varargin{k},'a_merra')]);     
   elseif strcmp(varargin{k},'a_merra2')==1
       % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from WRF
       aps_flag=30;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_merra2')]); 
   elseif strcmp(varargin{k},'a_merra-h')==1
       % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from WRF
       aps_flag=31;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_merra-h')]); 
   elseif strcmp(varargin{k},'a_merra2-h')==1
       % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from WRF
       aps_flag=32;   
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_merra2-h')]); 
  elseif strcmp(varargin{k},'a_merra-w')==1
       % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from WRF
       aps_flag=33;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_merra-w')]); 
   elseif strcmp(varargin{k},'a_merra2-w')==1
       % aps topo correlated modis recalibrated (non-interpolated) plus a hydrostatic component from WRF
       aps_flag=34; 
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_merra2-w')]);   
   elseif strcmp(varargin{k},'a_gacos')==1
       % aps topo correlated GACOS  
       aps_flag=35; 
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_gacos')]);    
  elseif strcmp(varargin{k},'a_narr')==1
       % aps topo correlated NARR  
       aps_flag=36; 
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_narr')]);  
  elseif strcmp(varargin{k},'a_narr-h')==1
       % aps topo correlated NARR  
       aps_flag=37; 
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_narr-h')]); 
  elseif strcmp(varargin{k},'a_narr-w')==1
       % aps topo correlated NARR  
       aps_flag=38; 
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_narr-w')]);        
  elseif strcmp(varargin{k},'a_era5')==1 || strcmp(varargin{k},'a_e5')==1
       % aps topo correlated ERA-I  correction
       aps_flag=39;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_era5') ]);
   elseif strcmp(varargin{k},'a_era5-h')==1 || strcmp(varargin{k},'a_e5h')==1
       % aps hydrostatic ERA-I correction
       aps_flag=40;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_era5-h') ]);
   elseif strcmp(varargin{k},'a_era5-w')==1 ||  strcmp(varargin{k},'a_e5w')==1
       % aps topo correlated ERA-I correction
       aps_flag=41;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_era5-w') ]);

  elseif strcmp(varargin{k},'i_as')==1
       % ionopsheric correlated atmopshere - azimuth shift method
       iono_flag=1;
       prmsrch=logical([prmsrch strcmp(varargin{k},'i_as')]);   
 
       
 elseif size(varargin{k},2)>=3  && strcmp(varargin{k}(1:3),'ifg')==1
           if size(varargin{k},2)==3
              ifg_number=[];
           else
              ifg_number = str2num(varargin{k}(5:end));
           end
           % display all bands for 'a_p' option for the ith interferogram
           aps_band_flag=1;
           prmsrch=logical([prmsrch strcmp(varargin{k}(1:3),'ifg') ]);  
   elseif size(varargin{k},2)>=3  && strcmp(varargin{k}(1:3),'ext')==1
           % Plot of external data
           if size(varargin{k},2)==3
               ext_data_path=[];
           else
              ext_data_path = varargin{k}(5:end);
           end
           ext_data_flag = 1;
           prmsrch=logical([prmsrch strcmp(varargin{k}(1:3),'ext') ]);          
   else
       prmsrch=logical([prmsrch 0]);
   end     
   
   
   
end

if sum(prmsrch) >= 1
  prmix=find(prmsrch);
  stdargin=stdargin-sum(prmsrch);  % drop one optional prm form the list of arguments
  varargin(prmsrch)=[]; % drop 'ts' from the list.
end


for k = 1:length(varargin(:))
    eval([arglist{k},'=varargin{',num2str(k),'};']);     
    %val=eval(['varargin{',num2str(k),'};']);         % debug
    %disp([arglist{k},'=varargin{',num2str(k),'};']);  % debug
    %disp([arglist{k}, num2cell(val)]);                % debug
end

% fix in case the specied APS option does not exist it will be mapped in
% plot_flag. However that needs to be a numeric value not a string. When
% this occurs give error as this plotting APS type does not exist
if ~isnumeric(plot_flag)
    error('This plotting option does not exist')
end

if plot_flag > 1 && ts_flag==1
    disp('TS plot is possible with backgrounds options 0 or 1')
    plot_flag=1
end

if plot_flag < 0 && ts_flag==1
    error('No time series plotting is possible with backgrounds option -1')
end
% checking if a valid option is slected for the aps_bands potting
if aps_band_flag==1 && aps_flag~=2 && aps_flag~=11
   error('myApp:argChk', ['For the aps display of all spatial bands only the "a_p" and "a_pk" flag can be selected. \n'])  
end
if aps_band_flag==1 && isempty(ifg_number)
   error('myApp:argChk', ['Specify the ith interferogram for the spatial bands option as "ifg i". \n'])  
end
% checking when selected to plot external data, the path exists
% if ext_data_flag==1 && (exist(ext_data_path,'dir')~=7 || exist(ext_data_path,'file')~=2)
%    error('myApp:argChk', ['External datapath folder/file does not exist. \n'])  
% end


%%% debug
%disp('summary: ')
%stdargin
%if Noptargin > 1
% plot_flag  % assigned later in the script that is called.
%end

%if ts_flag==0
%    disp('no ts_flag')
%end


