function []=ps_info()
%PS_INFO display information on the dataset
%
%   See also SB_INFO
%
%   Andy Hooper, March 2007
%
%   ======================================================================
%   01/2010 AH: Give single master info for small baselines too
%   ======================================================================


load psver
psname=['ps',num2str(psver)];

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

for i=1:size(ps.day,1)
    aa=[datestr(ps.day(i))];
    fprintf('%3s  %s %5s m\n',num2str(i),aa,num2str(round(bperp(i))));
end
fprintf('Number of stable-phase pixels: %d\n',ps.n_ps);

