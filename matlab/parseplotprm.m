% parse parameters for PS_PLOT function
%
%
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   11/2010 MA: initial support for 'ts' option 
%   ======================================================================

% list of arguments for ps_plot excluding value type
arglist={'plot_flag','lims','ref_ifg','ifg_list','n_x','cbar_flag',...
    'textsize','textcolor','lon_rg','lat_rg'};  % excluding value_type

%varargin  % debug
%stdargin=nargin ; 
%fprintf('Number of inputs = %d\n', nargin)

Noptargin = length(varargin(:));  % treat extra options as varargin, 
               % for ps_plotTS('v-d',4,'ts') is 3
ts_flag=0;

% search for char parameter like 'ts'
prmsrch=[];

for k=1:Noptargin,
   %prmsrch=[prmsrch ischar(varargin{k})];
   prmsrch=logical([prmsrch strcmp(varargin{k},'ts') ]);
end

if sum(prmsrch) >= 1
  prmix=find(prmsrch);
  ts_flag=1;
  stdargin=stdargin-1;  % drop one optional prm form the list of arguments
  varargin(prmsrch)=[]; % drop 'ts' from the list.
end
%length(varargin(:))

for k = 1:length(varargin(:))
    eval([arglist{k},'=varargin{',num2str(k),'};']);     
    %val=eval(['varargin{',num2str(k),'};']);         % debug
    %disp([arglist{k},'=varargin{',num2str(k),'};']);  % debug
    %disp([arglist{k}, num2cell(val)]);                % debug
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

