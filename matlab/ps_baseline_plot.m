function []=ps_baseline_plot
%PS_BASELINE_PLOT plot single master baselines
%
%   Andy Hooper, March 2011
%
%   ======================================================================
%   ======================================================================


load psver
psname=['ps',num2str(psver)];
ps=load(psname);

small_baseline_flag=getparm('small_baseline_flag');
if ~strcmpi(small_baseline_flag,'n')
   error('You are not in a PS directory')
end

m=ps.master_ix;


clf

for i=1:length(ps.day)
    l=line([ps.day(i),ps.day(m)],[ps.bperp(i),ps.bperp(m)]);
    set(l,'color',[0 1 0],'linewidth',2)
end
hold on
p=plot(ps.day,ps.bperp,'ro');
set(p,'markersize',12,'linewidth',2)
hold off
datetick('x',12)
ylabel('B_{perp}')


