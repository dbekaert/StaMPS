function [] = profiler(data,lonlat,point1,point2,bins,R,units)



% By David Bekaert - December 2010

% give as input the data, the lonlat (2 colums), the profile points [lon lat], the number of
% bins and the radius [m] of the tube used to project points
% by default a 15 bins and a radius of 100m is used.
% for debug figures set debugflag to 1

debugflag=0;


if nargin<7
    units = 'rad';
end
if nargin<6
    R = 100;            % default value of R in m
end
if nargin<5
    bins = 15;          % default number of bins
end
if nargin<4
    error('Specify more input arguments \n');    
end
fontsize = 15;          % label fontsize


%%% PLOTTING
profile_fig = figure('name','Cross-sectional plot');
axes('OuterPosition',[0  0  0.7  1])
scatter(lonlat(:,1),lonlat(:,2),3,data,'filled')
hold on
plot([point1(1,1) point2(1,1)]',[point1(1,2) point2(1,2)]','k-*','LineWidth',1.5)
cc =colorbar;
xlabel(cc,units,'fontsize',fontsize)
set(gca,'fontsize',fontsize)
xlabel('Longitude','fontsize',fontsize)
ylabel('Latitude','fontsize',fontsize)
%%% PLOTTING



% crossection computation
% transform to a local reference frame, unit becomes [km]
lonlatXY = llh2local(lonlat',point1);	% transformation of the PS pixels
P1 = llh2local(point1',point1);                 % transformation of the line points 
P2 = llh2local(point2',point1);

% rotate the local from such that line becomes horizontal
alpha = atan((P2(2,1)-P1(2,1))/(P2(1,1)-P1(1,1)));		% angle [rad]
rot = [cos(alpha)    -sin(alpha)
       sin(alpha)    cos(alpha)];
XY_new = rot\lonlatXY;                  % X is the first row, Y the second row
P1_new = rot\P1;		
P2_new = rot\P2;
clear rot alpha  P1 P2 lonlatXY



if debugflag==1
    figure('name','Debug_plot [1]');
    scatter(XY_new(1,:),XY_new(2,:),3,data,'filled')
    hold on
    scatter([P1_new(1) P2_new(1)],[P1_new(2) P2_new(2)],'ro')
end




% Search for all the PS within distance R
ix_data = find((abs(XY_new(2,:))<=R/1000));  

% updating vectors
XY_new = XY_new(:,ix_data);
data_new = data(ix_data);

if debugflag==1
    figure('name','Debug_plot [2]');
    scatter(XY_new(1,:),XY_new(2,:),3,data_new,'filled')
    hold on
    scatter([P1_new(1) P2_new(1)],[P1_new(2) P2_new(2)],'ro')
end


% output to the screen
PS_used = size(ix_data,2);
if (isempty(PS_used)==1 || PS_used==0)
    error('No PS are found within the tube crossection. \n')
else
    fprintf([num2str(PS_used),' PS used to compute projection on crossection \n'])
end

% threshold for minimum PS inside a bin
ps_min = floor(0.5*(PS_used/bins));

% binning of the results
binsize = (max(XY_new(1,:))-min(XY_new(1,:)))/bins;
data_cross_binned = zeros([1 bins]);
bins_xy = zeros([1 bins]);
clear ix

for k=1:bins
    bl = (k-1)*binsize+min(XY_new(1,:));                        % lower bound
    bu = (k)*binsize+min(XY_new(1,:));                          % upper bound
    ix = find(bl<=XY_new(1,:)  & XY_new(1,:)<bu);
    % Only showing bins with a minimum of PS contained    
    if size(ix,2)>=ps_min
        bins_xy(1,k) = bl+binsize/2;                            % position center of bin
        data_binned_mean(1,k) = mean(data_new(ix,1));     % take mean value
      
    else
        bins_xy(1,k) = NaN;
        data_binned_mean(1,k) = NaN;
        
    end    
    clear ix bl bu
end
ix = find(isnan(bins_xy)==1);
if isempty(ix)==0
    bins_xy(ix)=[];
    data_binned_mean(ix)=[];
end



%%% PLOTTING
figure(profile_fig)
axes('OuterPosition',[0.65  0  0.3  1])
plot(data_new,XY_new(1,:),'.','color',[0.9 0.9 0.9])
ylabel('Distance along crossection [km]','fontsize',fontsize)
hold on
plot(data_binned_mean,bins_xy,'k.')
xlabel(units,'fontsize',fontsize)
set(gca,'fontsize',fontsize)
%%% PLOTTING


clear all


