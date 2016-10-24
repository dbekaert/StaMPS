% ts_plot for PS_PLOT function
%
%
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   11/2010 MMC & MA: plotting time series using 'ts' option 
%   03/2011 AH convert radius to m, remove line fitting and add mean
%   08/2016 AH Replace break with return
%   ======================================================================

% Place button for new TS plots
% see ps_plot

mEditBox=findobj('Background',[1 1 1]); % get the handle of editable textbox
%radiusfactor=str2num(get(mEditBox,'String')) % get radius factor from editbox
radiusfactor = str2num(char(get(mEditBox,'String')));					% added by david this is also not a single value but a vector
radiusfactor = radiusfactor(1);								% select only the first value, i checked it for a couple of values and radius is correct

%if radiusfactor > 3600
%    disp('radius factor should be <= 3600')
%    return
%end

momfig_name=['(', get(gcf,'name'), ')']; % inherent fig_name from the velocity plot

% LOAD TS mat file ~ ps_plot_ts_v-d.mat
if ~exist('ph_mm','var')
  %fid = fopen(savetxt);
  fid = fopen('ps_plot_ts_matname.txt');
  tsmat=textscan(fid,'%s'); % get mat filename to load parameters
  fclose(fid);
  clear fid
  eval(['load ' tsmat{1}{1}])% load saved matrix
end
% end of load  
  
%  
% GET USER INPUT from MOUSE CLICK
%
 disp('Please select a point on the figure to plot time series (TS)')
 [lon0,lat0] = ginput(1); % turn on when final
 disp(['Selected point coordinates (lon,lat):' num2str(lon0),', ', num2str(lat0) ])

% MAKE A CIRCLE AROUND SELECTED POINT (lon0,lat0)
t = linspace(0,2*pi,114); % 114 pts
h=lon0; k=lat0;  % center cn for the circle
r=radiusfactor;
                   % change radius if you want to include more points
x = r*cos(t);
y = r*sin(t);

% SELECT POINTS based on LONLAT
xy=llh2local(lonlat',[lon0,lat0])'*1000;
in = inpolygon(xy(:,1),xy(:,2),x,y);
n_pts_near=sum(in);  % how many ps found

if sum(in) == 0
    disp(['No points found in radius of ', num2str(r) ] )
    disp('Please make new selection increasing radius factor...')
    return     % no pts selecte
end

% PLOT selection and points
figure
    set(gcf,'name',['Found #pt(s): ', num2str(n_pts_near),...
        '  in radius: ', num2str(r), ' m. ', momfig_name ])
    circlell=local2llh([x;y]/1000,[lon0,lat0]);
    plot(circlell(1,:),circlell(2,:));
    hold on
    plot(lon0,lat0,'*')

    lon2=lonlat(in,1);
    lat2=lonlat(in,2);
    plot(lon2,lat2,'dr')               
    axis image; hold off
    
% if ps>1 than average
% [dist,az] = distance(lat0,lon0,lat2,lon2); % or use llh2local

% plot closest point or avg of multiple points.


% phases
% ph_mm holds corrected phase
% ph_all holds radians to meters
% v_all=-m(2,:)'; % ?
  
% PLOT TS for given point(s)
ts=ph_mm(in,:);
G=[ones(size(day)),day-master_day] ; % [ 1  a ] --> b + ax

offset=pi*1000*lambda/(4*pi);

ts=nanmean(ts,1);
x_hat=G\double(ts');
ts_hat=G*x_hat;
%tsup_hat=ts_hat+offset;
%tslo_hat=ts_hat-offset;


%%% Typical TS plot - no auxilary subplots
% figure
%     set(gcf,'name',[ ' Times series plot for #point(s): ',...
%         num2str(n_pts_near), ' ', momfig_name])
%     plot(day,ts,'--*'); hold on
%     plot(day,ts_hat,'-*r');
%     plot(day,tsup_hat,'-.g');
%     plot(day,tslo_hat,'-.g');
%     hold off 
%     grid on
%     ylabel('mm');
%     xlabel('Time [mmmyy]')
%     %datetick('x','mmmyy')  % keepticksdoc 
%     %set(gca, 'XTick',day);
%     set(gca, 'XTickLabel', datestr(day,'mmmyy'));

%%% Enhanced TS plot
 figure % main figure
   orient landscape
   % Bperp
%    subplot(10,1,1)     % Bperp
%      bperp(find(bperp==0))=[]; % drop master.
%      bar(bperp)
%      ylabel('Bperp [m]')
%    grid on
   % TS
%    subplot(10,10,[11 87]) % subplot(10,1,2:9)
       set(gcf,'name',[ ' Times series plot for #point(s): ',...
        num2str(n_pts_near), ' ', momfig_name])
    %h1=plot(day,ts,'--*'); hold on
    %plot(day,ts_hat,'-k','LineSmoothing','on'); % mess up ticks
    plot(day,ts_hat,'-r'); % mess up ticks
    hold on
    h2=plot(day,ts,'ob');
    set(h2,'linewidth',2)
    %h3=plot(day,tsup_hat,'-.g');
    %h4=plot(day,tslo_hat,'-.g');
    set(gca,'fontsize',20)
    hold off
    grid on
    ylabel('LOS (mm)');
    %datetick('x','mmmyy')  % keepticks or see below
   
    ts_years=unique(datenum(datestr(day,'yyyy'),'yyyy'));
    ts_dates=[ts_years(1):365.25:ceil(ts_years(end)+365.25)];
    %ts_dates=[day(1):365:day(end)+365];
    set(gca, 'XTick',ts_dates);
    set(gca, 'XTickLabel', datestr(ts_dates,'yyyy'));
    % annotate velocit slope in ts plot  
     
   % IFG Dates - excluding master 
%    subplot(10,10,[18 90]) % subplot for rectangle
%       putdates(0.05,1,datestr(day,'yyyy-mm-dd'),0.035,9) % putdates(xstart, ystart, labels, labeloffset, fontsize)
   %subplot(10,1,10)
   %bar(Bdop)
    
       
%EOF
