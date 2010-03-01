%BATCHJOB initial script invoked for all background jobs
%
%   Andy Hooper, Sep 2001

[p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10),p(11),p(12),therest]=textread('matbgparms.txt','%s%s%s%s%s%s%s%s%s%s%s%s%[^\n]');

%if funcname(end-1)=='.' 
%  funcname=funcname(1:end-2);
%end

parmarray{1}=char(p{1});

j=1;
for i=2:length(p)
    switch isempty(p{i})
    case 0
        j=j+1;
        pnum=str2num(p{i});
        switch isempty(pnum)
        case 1
             parmarray{j}=p{i};
        otherwise
             parmarray{j}=pnum;
        end
    end
end

parmarray
!hostname

%setpath

tic
feval(parmarray{:});
toc
