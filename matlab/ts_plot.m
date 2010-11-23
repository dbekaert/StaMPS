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
% mButton=uicontrol('Style', 'pushbutton', 'Callback', 'ts_plot',...
%     'String','new TS plot', 'Position', [240 30 90 30] ); % if exist don't create
% 
% mButtonposition=get(mButton,'Position')
% 
% mTextBox=uicontrol('Style', 'edit',...
%     'String','1', 'Position', [280 30 90 30] );

mEditBox=findobj('Background',[1 1 1]); % get the handle of editable textbox
%radiusfactor=str2num(get(mEditBox,'String')) % get radius factor from editbox
radiusfactor=str2num(char(get(mEditBox,'String')));					% added by david this is also not a single value but a vector
radiusfactor = radiusfactor(1);								% select only the first value, i checked it for a couple of values and radius is correct

if radiusfactor > 100
    disp('radius factor should be <= 100')
    break
end

% if ph_uw does not exist load otherwise don't load TODO
momfig_name=['(', get(gcf,'name'), ')']; % inherent fig_name from the velocity plot

% clear all % clean up
  %fid = fopen(savetxt);
  fid = fopen('ps_plot_ts_matname.txt');
  tsmat=textscan(fid,'%s'); % get mat filename to load parameters
  fclose(fid);
  clear fid
  eval(['load ' tsmat{1}{1}])% load saved matrix

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
%r=1/3600*2;        % 1 arcsec: 1/3600 radius --> 0.000277777777777778 ~ 30 m at the equator    
r=1/3600*radiusfactor;
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
    lat2=lonlat(in,2);
    plot(lon2,lat2,'dr')               
    axis image; hold off
    
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
x_hat=G\double(ts');

offset=pi*1000*lambda/(4*pi);
x1_hat=[x_hat(1)+offset; x_hat(2)];
x2_hat=[x_hat(1)-offset; x_hat(2)];

ts_hat=G*x_hat;
tsup_hat=G*x1_hat;
tslo_hat=G*x2_hat;
% whos ts G x_hat ts_hat

% % typical TS plot
% figure
%     set(gcf,'name',[ ' Times series plot for #point(s): ',...
%         num2str(n_pts_near), ' ', momfig_name])
%     plot(day./10^5,ts,'--*'); hold on
%     plot(day./10^5,ts_hat,'-*r','LineSmoothing','on');
%     plot(day./10^5,tsup_hat,'-.g');
%     plot(day./10^5,tslo_hat,'-.g');
%     hold off  % this is crucial for datetick 
%     grid on
%     ylabel('mm');
%     xlabel('Time [mmmyy]')
%     format long g
%     %datetick('x','mmmyy')  % keepticksdoc 
%     %set(gca, 'XTick',day./10^5);
%     set(gca, 'XTickLabel', datestr(day,'mmmyy'));

% enhanced TS plot
 figure % main figure
   %orient landscape
   subplot(10,1,1)     % Bperp
     bperp(find(bperp==0))=[]; % drop master.
     bar(bperp)
     ylabel('Bperp [m]')
   grid on
   
   %title(['Average coherence of ' num2str( n_ps ) ' PSs with coherence window ' num2str( win_size(1) ) 'x' num2str( win_size(2) ) ' arranged in BDop']);
   subplot(10,10,[11 87]) % subplot(10,1,2:9)
       set(gcf,'name',[ ' Times series plot for #point(s): ',...
        num2str(n_pts_near), ' ', momfig_name])
    h1=plot(day,ts,'--*'); hold on
    %plot(day,ts_hat,'-*r','LineSmoothing','on'); % mess up ticks
    h2=plot(day,ts_hat,'-*r');
    h3=plot(day,tsup_hat,'-.g');
    h4=plot(day,tslo_hat,'-.g');
    hold off
    grid on
    ylabel('mm');
    xlabel('Time [mmmyy]')
    %format long g
    datetick('x','mmmyy')  % keepticksdoc 

    %set(gca, 'XTick',day./10^5);
    %set(gca, 'XTickLabel', datestr(day,'mmmyy'));
    
   subplot(10,10,[18 90]) % subplot for rectangle
      putdates(0.05,1,datestr(day,'yyyy-mm-dd'),0.035,9) % putdates(xstart, ystart, labels, labeloffset, fontsize)
   %subplot(10,1,10)
   %bar(Bperp(ind_Bdfdc))
    

    
% annotate slope    
    
%EOF
