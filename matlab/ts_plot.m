% ts_plot for PS_PLOT function
%
%
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   11/2010 MMC & MA: plotting time series using 'ts' option 
%   ======================================================================

% Place button for new TS plots
uicontrol('Style', 'pushbutton', 'Callback', 'ts_plot',...
    'String','new TS plot', 'Position', [240 30 90 30] ); % if exist don't create

% if ph_uw does not exist load otherwise don't load TODO
momfig_name=['(', get(gcf,'name'), ')']; % inherent fig_name from the velocity plot

whos
% clear all % clean up
  %fid = fopen(savetxt);
  fid = fopen('ps_plot_ts_matname.txt')
  tsmat=textscan(fid,'%s'); % get mat filename to load parameters
  fclose(fid);
  clear fid
  eval(['load ' tsmat{1}{1}])% load saved matrix

% end of load  
  
%  
% GET USER INPUT from MOUSE CLICK
%
 disp('Please select a point on the figure to plot time series (TS)')
 [lon0,lat0] = ginput(1) % turn on when final
 disp(['Selected point coordinates (lon,lat):' num2str(lon0),', ', num2str(lat0) ])

% MAKE A CIRCLE AROUND SELECTED POINT (lon0,lat0)
t = linspace(0,2*pi,114); % 114 pts
h=lon0; k=lat0;  % center cn for the circle
r=1/3600*2;        % 1 arcsec: 1/3600 radius --> 0.000277777777777778 ~ 30 m at the equator    
                   % change radius if you want to include more points
x = r*cos(t)+h;
y = r*sin(t)+k;
% for inpolygon repeat first cn at the end
xv = [x x(1)]; yv = [y y(1)];

% SELECT POINTS based on LONLAT
in = inpolygon(lonlat(:,1),lonlat(:,2),xv,yv);
n_pts_near=sum(in);  % how many ps found

if sum(in) == 0
    disp(['No points found in radius of ', num2str(r) ] )
    disp('Please make new selection...')
    break     % no pts selecte
end


figure
    set(gcf,'name',['Found #pt(s): ', num2str(n_pts_near),...
        '  in radius: ', num2str(r), ' deg. ', momfig_name ])
    plot(x,y);
    hold on
    plot(lon0,lat0,'*')

    lon2=lonlat(in,1);
    lat2=lonlat(in,2)
    plot(lon2,lat2,'dr')               
    axis image;
    
% if ps>1 than average
[dist,az] = distance(lat0,lon0,lat2,lon2); % or use llh2local


% plot closest point or avg of multiple points.


% phases
% ph_uw holds corrected phase
% ph_all holds radians to meters
%v_all=-m(2,:)'; % ?
  
% PLOT TS for given point(s)
ts=-ph_uw(in,:)*lambda*1000/(4*pi);
G=[ones(size(day)),day-master_day] ; % [ 1  a ] --> b + ax
x_hat=G\double(ts')

offset=pi*1000*lambda/(4*pi);
x1_hat=[x_hat(1)+offset; x_hat(2)];
x2_hat=[x_hat(1)-offset; x_hat(2)];

ts_hat=G*x_hat;
tsup_hat=G*x1_hat;
tslo_hat=G*x2_hat;
whos ts G x_hat ts_hat

figure
    set(gcf,'name',[ ' Times series plot for #point(s): ',...
        num2str(n_pts_near), ' ', momfig_name])
    plot(day,ts,'--*'); hold on
    plot(day,ts_hat,'-*r','LineSmoothing','on');
    plot(day,tsup_hat,'-.g');
    plot(day,tslo_hat,'-.g');
    grid on
    ylabel('mm');
    datetick('x','mmmyy')  

% annotate slope    
    
%EOF