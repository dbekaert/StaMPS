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
%   ======================================================================

load psver
psname=['ps',num2str(psver)];
ifgstdname=['ifgstd',num2str(psver)];

small_baselines_flag=getparm('small_baselines_flag');

ps=load(psname);

if isfield(ps,'ifgday')
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

if ~exist([ifgstdname,'.mat'],'file') & ~strcmpi(small_baselines_flag,'y')
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
    % fprintf('%3s  %s %5s m\n',num2str(i),aa,num2str(round(bperp(i))));
    fprintf('%3s  %s %5s m   %.3f deg\n',num2str(i),aa,num2str(round(bperp(i))), ifg_std(i));
end
fprintf('Number of stable-phase pixels: %d\n',ps.n_ps);

