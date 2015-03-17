function []=ps_baseline_plot(linecolor)
%PS_BASELINE_PLOT plot single master baselines
%
%   Andy Hooper, March 2011
%
%   ======================================================================
%   ======================================================================
if nargin<1
    linecolor=[1 0 0];
end

load psver
psname=['ps',num2str(psver)];
ps=load(psname);

small_baseline_flag=getparm('small_baseline_flag');
if ~strcmpi(small_baseline_flag,'n')
   error('You are not in a PS directory')
end

m=ps.master_ix;


%clf

for i=1:length(ps.day)
    l=line([ps.day(i),ps.day(m)],[ps.bperp(i),ps.bperp(m)]);
    set(l,'color',linecolor,'linewidth',2)
end
hold on
set(gca,'fontsize',14)
p=plot(ps.day,ps.bperp,'k+');
set(p,'markersize',8,'linewidth',1)
hold off
datetick('x',10)
xlabel('Acquisition Date')
ylabel('Perpendicular Baseline (m)')


