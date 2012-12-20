function []=sb_info()
%SB_INFO display small baseline information on the dataset
%
%   Andy Hooper, January 2010
% Modifications:
% 12/2012   DB:     Add std information as output

load psver
psname=['ps',num2str(psver)];
ps=load(psname);

% loading the std information
ifgstdname=['ifgstd',num2str(psver)];
if exist([ifgstdname,'.mat'],'file')
    stdin=load(ifgstdname);
    ifg_std=stdin.ifg_std;
    clear stdin
else
    ifg_std=zeros(ps.n_ifg,1);
end

if isfield(ps,'ifgday')
    for i=1:ps.n_ifg
        aa=[datestr(ps.ifgday(i,1)),' to ',datestr(ps.ifgday(i,2))];
        fprintf('%3s  %s %5s m %.3f deg\n',num2str(i),aa,num2str(round(ps.bperp(i))),ifg_std(i));
    end
    fprintf('Number of stable-phase pixels: %d\n',ps.n_ps);
else
    fprintf('This is not a small baseline directory\n');
end

