function []=ps_load_info()
%PS_LOAD_INFO Initial load of PS info files into matlab workspace

fprintf('Loading info into matlab...\n')

bperpname=['bperp.1.in'];           % in meters 1 line per ifg
dayname=['day.1.in'];               % YYYYMMDD, 1 line per ifg
masterdayname=['master_day.1.in'];  % YYYYMMDD
headingname=['heading.1.in'];       % satellite heading
lambdaname=['lambda.1.in'];         % wavelength
widthname=['width.txt'];            % width of interferograms
lenname=['len.txt'];                % length of interferograms

psver=1;

if ~exist(dayname,'file')
    dayname= ['../',dayname];
end
day=load(dayname);
year=floor(day/10000);
month=floor((day-year*10000)/100);
monthday=day-year*10000-month*100;
slave_day=datenum(year,month,monthday);
[slave_day,day_ix]=sort(slave_day);

if ~exist(masterdayname,'file')
    masterdayname= ['../',masterdayname];
end
master_day=load(masterdayname);
master_day_yyyymmdd=master_day;
year=floor(master_day/10000);
month=floor((master_day-year*10000)/100);
monthday=master_day-year*10000-month*100;
master_day=datenum(year,month,monthday);

master_ix=sum(slave_day<master_day)+1;
day=[slave_day(1:master_ix-1);master_day;slave_day(master_ix:end)]; % insert master day
n_ifg=length(day);

if ~exist(bperpname,'file')
    bperpname= ['../',bperpname];
end
bperp=load(bperpname); 
bperp=bperp(day_ix);
bperp=[bperp(1:master_ix-1);0;bperp(master_ix:end)]; % insert master-master bperp (zero)
n_ifg=size(bperp,1);

 
if ~exist(headingname,'file')
    headingname= ['../',headingname];
end
heading=load(headingname);
if isempty(heading)
    error('heading.1.in is empty')
end
setparm('heading',heading,1);

if ~exist(lambdaname,'file')
    lambdaname= ['../',lambdaname];
end
lambda=load(lambdaname);
setparm('lambda',lambda,1);

n_ps=0;

savename=['ps',num2str(psver)];

save(savename,'n_ifg','bperp','day','master_day','master_ix','day_ix','n_ps');
save psver psver
