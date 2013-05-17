function []=ps_info()
%PS_INFO display information on the dataset
%
%   See also SB_INFO
%
%   Andy Hooper, March 2007
%
%   ======================================================================
%   01/2010 AH: Give single master info for small baselines too
%   10/2010 MA: added noise estimates 
%   06/2011 AH: remove noise estimates for small baselines
%   10/2011 AH: check for small baselines before inverting for bperp
%   04/2013 DB: Fixed type for small_baseline_flag
%   ======================================================================

load psver
psname=['ps',num2str(psver)];
ifgstdname=['ifgstd',num2str(psver)];
rcname=['rc',num2str(psver)];
pmname=['pm',num2str(psver)];

small_baseline_flag=getparm('small_baseline_flag');

ps=load(psname);
if strcmpi(small_baseline_flag,'y')
    G=zeros(ps.n_ifg,ps.n_image);
    for i=1:ps.n_ifg
         G(i,ps.ifgday_ix(i,1))=-1;
         G(i,ps.ifgday_ix(i,2))=1;
    end
    G=G(:,[1:ps.master_ix-1,ps.master_ix+1:end]);
    bperp=[G\double(ps.bperp)];
    bperp=[bperp(1:ps.master_ix-1);0;bperp(ps.master_ix:end)];
else
    bperp=ps.bperp;
end

if ~exist([ifgstdname,'.mat'],'file') & ~strcmpi(small_baseline_flag,'y') & exist([rcname,'.mat'],'file') & exist([pmname,'.mat'],'file')
    ps_calc_ifg_std
end

if exist([ifgstdname,'.mat'],'file')
    stdin=load(ifgstdname);
    ifg_std=stdin.ifg_std;
    clear stdin
else
    ifg_std=zeros(size(ps.day,1),1);
end

for i=1:size(ps.day,1)
    aa=[datestr(ps.day(i))];
    if ifg_std(i)==0
      fprintf('%3s  %s %5s m\n',num2str(i),aa,num2str(round(bperp(i))));
    else
      fprintf('%3s  %s %5s m   %.3f deg\n',num2str(i),aa,num2str(round(bperp(i))), ifg_std(i));
    end 
end
fprintf('Number of stable-phase pixels: %d\n',ps.n_ps);

