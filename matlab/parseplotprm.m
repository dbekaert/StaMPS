%   parse parameters for PS_PLOT function
%
%   Andy Hooper, June 2006
%
%   modifications:
%   11/2010 MA: initial support for 'ts' option 
%   04/2013 DB: include support of topocorrelated aps 'a_m' for merris, 
%               'a_l' for linear, and 'a_p' for powerlaw option
%   04/2013 DB: include bandfilter option 'ifg i'

% list of arguments for ps_plot excluding value type
arglist={'plot_flag','lims','ref_ifg','ifg_list','n_x','cbar_flag',...
    'textsize','textcolor','lon_rg','lat_rg'};  % excluding value_type

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
plot_flag=1; 

% search for char parameter like 'ts','a_l','a_m','a_p'
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
       
       
   elseif size(varargin{k},2)>3  && strcmp(varargin{k}(1:3),'ifg')==1
           % display all bands for 'a_p' option for the ith interferogram
           aps_band_flag=1;
           ifg_number = str2num(varargin{k}(5:end));
           prmsrch=logical([prmsrch strcmp(varargin{k}(1:3),'ifg') ]);   
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



%%% debug
%disp('summary: ')
%stdargin
%if Noptargin > 1
% plot_flag  % assigned later in the script that is called.
%end

%if ts_flag==0
%    disp('no ts_flag')
%end

