function []=sb_info()
%SB_INFO display small baseline information on the dataset
%
%   Andy Hooper, January 2010

load psver
psname=['ps',num2str(psver)];

ps=load(psname);

if isfield(ps,'ifgday')
    for i=1:ps.n_ifg
        aa=[datestr(ps.ifgday(i,1)),' to ',datestr(ps.ifgday(i,2))];
        fprintf('%3s  %s %5s m\n',num2str(i),aa,num2str(round(ps.bperp(i))));
    end
    fprintf('Number of stable-phase pixels: %d\n',ps.n_ps);
else
    fprintf('This is not a small baseline directory\n');
end

