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
   elseif strcmp(varargin{k},'a_l')==1
       % aps topo correlated linear correction
       aps_flag=1;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_l') ]);
   elseif strcmp(varargin{k},'a_p')==1
       % aps topo correlated powerlaw correction
       aps_flag=2;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_p') ]);
   elseif strcmp(varargin{k},'a_m')==1
       % aps topo correlated meris correction
       aps_flag=3;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_m') ]);
   elseif strcmp(varargin{k},'a_e')==1
       % aps topo correlated ERA-I  correction
       aps_flag=4;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_e') ]);
   elseif strcmp(varargin{k},'a_eh')==1
       % aps hydrostatic ERA-I correction
       aps_flag=5;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_eh') ]);
   elseif strcmp(varargin{k},'a_ew')==1
       % aps topo correlated ERA-I correction
       aps_flag=6;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_ew') ]);
   elseif strcmp(varargin{k},'a_w')==1
       % aps topo correlated WRF correction
       aps_flag=7;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_w') ]);
   elseif strcmp(varargin{k},'a_wh')==1
       % aps hydrostatic WRF correction
       aps_flag=8;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_wh') ]);
   elseif strcmp(varargin{k},'a_ww')==1
       % aps topo correlated WRF correction
       aps_flag=9;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_ww') ]);
   elseif strcmp(varargin{k},'a_mi')==1
       % aps topo correlated  MERIS (non-interpolated)
       aps_flag=10;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_mi') ]);     
  elseif strcmp(varargin{k},'a_pk')==1
       % Spatial maps of the coefficient relating phase and tropo for power-law
       aps_flag=11;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_pk') ]);   
 elseif strcmp(varargin{k},'a_M')==1
       % aps topo correlated modis correction
       aps_flag=12;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_M') ]);       
 elseif strcmp(varargin{k},'a_MI')==1
       % aps topo correlated modis (non-interpolated)
       aps_flag=13;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_MI') ]);
       
  elseif strcmp(varargin{k},'a_m+a_eh')==1
       % aps topo correlated MERIS plus a hydrostatic component from ERA-I
       aps_flag=14;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_m+a_eh')]);     
  elseif strcmp(varargin{k},'a_mi+a_eh')==1
       % aps topo correlated MERIS (non-interpolated) plus a hydrostatic component from ERA-I
       aps_flag=15;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_mi+a_eh')]);     
  elseif strcmp(varargin{k},'a_M+a_eh')==1
       % aps topo correlated modis plus a hydrostatic component from ERA-I
       aps_flag=16;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_M+a_eh')]);     
  elseif strcmp(varargin{k},'a_MI+a_eh')==1
       % aps topo correlated modis (non-interpolated) plus a hydrostatic component from ERA-I
       aps_flag=17;
       prmsrch=logical([prmsrch strcmp(varargin{k},'a_MI+a_eh')]);    
       
       
       
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
    
if plot_flag > 1 && ts_flag==1
    disp('TS plot is possible with backgrounds options 0 or 1')
    plot_flag=1
end

if plot_flag < 0 && ts_flag==1
    disp('No time series plotting is possible with backgrounds option -1')
    break
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


