function []=sb_find_delaunay(n_consec)
%SB_FIND_consec find small baselines
%   SB_FIND_CONSEC find the n_consec ifg combinations
%
%
%   David Bekaert, April 2017
%
%   ======================================================================
%   ======================================================================



load psver
psname=['ps',num2str(psver)];
ps=load(psname);


if nargin<1
    n_consec = 3;
end

dates_ix= [];
for date=1:length(ps.day)
    dates_ix_new = repmat(date,n_consec,2);
    dates_ix_new(:,2) = dates_ix_new(:,2) + [1:n_consec]';
    dates_ix = [dates_ix ;dates_ix_new];
end
% remove the dates which are outside the network
ix_drop = find(dates_ix(:,2)>length(ps.day));
dates_ix(ix_drop,:)=[];

sbname='small_baselines.list';


x=dates_ix(:,1);
y=dates_ix(:,2);

clf

for i=1:length(x)
    l=line([ps.day(x(i)),ps.day(y(i))],[ps.bperp(x(i)),ps.bperp(y(i))]);
    text((ps.day(x(i))+ps.day(y(i)))/2,(ps.bperp(x(i))+ps.bperp(y(i)))/2,num2str(i));
    set(l,'color',[0 1 0],'linewidth',2)
end
hold on
p=plot(ps.day,ps.bperp,'ro');
set(p,'markersize',12,'linewidth',2)
hold off
datetick('x',12)
ylabel('B_{perp}')


h=gcf;
print(h,'-dpng','SB_plot.png')


fid=fopen(sbname,'w');
for i=1:length(x)
    fprintf(fid,'%s %s\n',datestr(ps.day(x(i)),'yyyymmdd'),datestr(ps.day(y(i)),'yyyymmdd'));
end
fclose(fid);

