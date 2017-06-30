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


if nargin <1
   ix=[];
end

fontsize = 25;

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

n_ifg=size(sb,1);
[yyyymmdd,I,J]=unique(sb);
ifg_ix=reshape(J,n_ifg,2);
x=ifg_ix(:,1);
y=ifg_ix(:,2);


day=str2num(datestr(ps.day,'yyyymmdd'));
[B,I]=intersect(day,yyyymmdd);

x=I(x);
y=I(y);

h_fig=figure('name','Baseline plot','position',[ 440   170   821   628]);
for i=1:length(x)
    l=line([ps.day(x(i)),ps.day(y(i))],[ps.bperp(x(i)),ps.bperp(y(i))]);
    text((ps.day(x(i))+ps.day(y(i)))/2,(ps.bperp(x(i))+ps.bperp(y(i)))/2,num2str(ix(i)),'fontsize',fontsize-7);
    set(l,'color',[0 1 0],'linewidth',2)
end
hold on
p=plot(ps.day,ps.bperp,'ro');
set(p,'markersize',12,'linewidth',2)
hold off
box on
set(gca,'fontsize',fontsize)
datetick('x',12)
ylabel('B_{perp}','fontsize',fontsize)
try
    print(h_fig,'-dpng','baseline_plot.png')
catch
end


