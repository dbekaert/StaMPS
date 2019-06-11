function []=plot_sb_baselines(ix)
%PLOT_SB_BASELINES plot the small baselines in small_baselines.list
%Optional an input argument ix can be specified containing a vector of the
%small baseline interferograms to keep which will be plotted in the baseline plot.
%
%   Andy Hooper, June 2007
%
%   ======================================================================
%   09/2010 AH: Add option to plot in MERGED directory
%   09/2010 AH: For SMALL_BASELINES/MERGED don't plot dropped ifgs 
%   12/2012 DB: Added meaning of ix to the syntax of the code
%   04/2013 DB: Command variable
%   03/2014 DB: Suppress command line output
%   01/2017 DB: Relax the name for SB such it could have a prefix or suffix
%   01/2017 DB: add fontsize
%   ======================================================================

% 05/2019 MG: Modified to plot colorize network based on interferogram standard
% deviation 


if nargin <1
   ix=[];
end

fontsize = 16;

currdir=pwd;
dirs=strread(currdir,'%s','delimiter','/');


if ~isempty(findstr(dirs{end},'SMALL_BASELINES'))
    
    if strcmp(dirs{end},'SMALL_BASELINES') 
        [a,b] = system(['\ls -d [1,2]* | sed ''' 's/_/ /''' ' > small_baselines.list']);
    else
        cd ../SMALL_BASELINES
        [a,b] = system(['\ls -d [1,2]* | sed ''' 's/_/ /''' ' > ' currdir filesep 'small_baselines.list']);
        cd(currdir)
    end
    load ../psver
    psname=['../ps',num2str(psver)];
    small_baseline_flag='y';
 
elseif strcmp(dirs{end},'MERGED') 
    cd ../SMALL_BASELINES
    [a,b] = system(['\ls -d [1,2]* | sed ''' 's/_/ /''' ' > ../MERGED/small_baselines.list']);
    cd ../MERGED
    load ../psver
    psname=['../ps',num2str(psver)];
    small_baseline_flag='y';
else
    load psver
    psname=['ps',num2str(psver)];
    small_baseline_flag='n';
end

sb=load('small_baselines.list');
n_ifg=size(sb,1);
if small_baseline_flag=='y' & isempty(ix) & exist('./parms.mat','file')
    drop_ifg_index=getparm('drop_ifg_index');
    if ~isempty(drop_ifg_index)
       ix=setdiff([1:n_ifg],drop_ifg_index);
    end
end

if ~isempty(ix)
    sb=sb(ix,:);
else 
    ix=1:size(sb,1);
end

ps=load(psname);

% loading the std information
    ifgstdname=['ifgstd',num2str(2)];
    if exist([ifgstdname,'.mat'],'file')
       stdin=load(ifgstdname);
       ifg_std=stdin.ifg_std;
       clear stdin
    else
    ifg_std=zeros(ps.n_ifg,1);
    end
    
    if ~isempty(ix)
    ifg_std = ifg_std(ix,:);
    else 
    ix=1:size(ifg_std,1);
    end

n_ifg=size(sb,1);
[yyyymmdd,I,J]=unique(sb);
ifg_ix=reshape(J,n_ifg,2);
x=ifg_ix(:,1);
y=ifg_ix(:,2);


day=str2num(datestr(ps.day,'yyyymmdd'));
[B,I]=intersect(day,yyyymmdd);

x=I(x);
y=I(y);

%create colormap
%colors
C1 = [1 0 0]; %red
C2 =[1 1 1];  %white
C3 = [0 0 1];  %blue

invert_flag = 0;

%invert colormap
if invert_flag ==1
  CM1=[linspace(C1(1),C2(1),100)',linspace(C1(2),C2(2),100)', linspace(C1(3),C2(3),100)'];
  CM2=[linspace(C2(1),C3(1),100)',linspace(C2(2),C3(2),100)', linspace(C2(3),C3(3),100)'];
  CM = [CM1;CM2];
elseif invert_flag ==0
  CM1=[linspace(C2(1),C1(1),100)',linspace(C2(2),C1(2),100)', linspace(C2(3),C1(3),100)'];
  CM2=[linspace(C3(1),C2(1),100)',linspace(C3(2),C2(2),100)', linspace(C3(3),C2(3),100)'];
  CM = [CM2;CM1];
else
end

%cm = colormap('parula'); % take your pick (doc colormap)
%cm = CM;
%cm = interp1(linspace(min(ifg_std),max(ifg_std),length(cm)),cm,ifg_std); % map color to y values
%cm = uint8(cm'*255); % need a 4xN uint8 array
%cm(4,:) = 255; % last column is transparency

h_fig=figure('name','Baseline plot','position',[ 440   170   821   628]);
for i=1:length(x)
    hold on
    x1 = [ps.day(x(i)),ps.day(y(i))]';
    y1 = [ps.bperp(x(i)),ps.bperp(y(i))]';
    xx = [x1 x1];
    yy = [y1 y1];
    zz = zeros(size(xx));
    cc = [ifg_std(i) ifg_std(i)];
    l=surf(xx,yy,zz,[cc;cc],'EdgeColor','interp','FaceColor','none');
    %l=plot([ps.day(x(i)),ps.day(y(i))],[ps.bperp(x(i)),ps.bperp(y(i))],'-');
    text((ps.day(x(i))+ps.day(y(i)))/2,(ps.bperp(x(i))+ps.bperp(y(i)))/2,num2str(ix(i)),'fontsize',fontsize-7);
    set(l,'linewidth',2)
    %set(l.Edge,'ColorBinding','interpolated','ColorData',cm)
    %drawnow i
end

if length(ix) ~= 1    
colormap(CM);
cb = colorbar;
ylabel(cb,'Ifg std [deg]');
else
end
shading flat
view(2)
p=plot(ps.day,ps.bperp,'ko');
set(p,'markersize',6,'linewidth',2)
hold off
box on
set(gca,'fontsize',fontsize)
datetick('x',12)
ylabel('B_{perp}','fontsize',fontsize)
try
    print(h_fig,'-dpng','baseline_plot.png')
catch
end


