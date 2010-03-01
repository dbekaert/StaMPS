function []=plot_all_ifgs(type,n_x,textsize)
% PLOT_ALL_IFGS plot all interferograms
%
%   Andy Hooper, June 2006
%
%   ======================================================
%   09/2006 AH:  small baselines added 
%   03/2007 AH:  variable text size added
%   04/2007 AH:  azimuth interferograms added
%   09/2008 AH:  generalise for other sensors
%   ======================================================


if nargin<1
    type=0;
end

if nargin<2
    n_x=0;
end

if nargin<3
    textsize=0;
end




heading=getparm('heading');
if isempty(heading)
   heading=-160;
end

if exist('looks.txt','file')
    looks=load('looks.txt');
else 
    looks=4;
end

if type==0
    ifgfile=['cint.minrefdem_',num2str(looks),'l.ras'];
else
    ifgfile=['cint.azint_',num2str(looks),'l.ras'];
end

%y=linspace(0.80,0.05,5);
%x=linspace(0.05,0.85,6);
%[imY,imX]=meshgrid(y,x);
%lx=0.14;
%ly=0.135;
i2=0;
ifgstack=[];
ifgname=[];
ifgname2=[];

aa=dir;
for i=1:size(aa,1)
    if aa(i).isdir==1 && strncmp(aa(i).name,'.',1)~=1
        dirname=aa(i).name;
        %if ~strncmp(aa(i).name,'1992',4) & ~strncmp(aa(i).name,'1993',4)
        cd(aa(i).name)
        if exist(ifgfile,'file')
            i2=i2+1;
            [ifg,cc]=imread(ifgfile);
            ifgstack(:,:,i2)=ifg;
            ifgname(i2)=str2num(dirname(1:8));
            if length(dirname)==17
                ifgname2(i2)=str2num(dirname(10:17));
            end
        end
        cd ..
        %end
    end
end

bb=dir('*crop.slc');
if ~isempty(bb)
   mastername=str2num(bb.name(1:8));
   master_ix=sum(mastername>ifgname)+1;
   ifgname=[ifgname(1:master_ix-1),mastername,ifgname(master_ix:end)];
   ifgstack=cat(3,ifgstack(:,:,1:master_ix-1),ones(size(ifgstack(:,:,1)))*113,ifgstack(:,:,master_ix:end));
else 
    master_ix=0;
end

[n_i,n_j,n_k]=size(ifgstack);
fig_ar=1.33;
ar=n_j/n_i/fig_ar;
if n_x==0
    n_y=ceil(sqrt(n_k*ar)); % number of plots in y direction
    n_x=ceil(n_k/n_y);
else
    n_y=ceil(n_x*ar) % number of plots in y direction
end
%n_y=ceil(sqrt(n_k*ar)); % number of plots in y direction
%n_x=ceil(n_k/n_y);
d_y=1/n_y;
d_x=d_y*ar;
if d_x>1/n_x
   d_x=1/n_x;
   d_y=d_x/ar;
end

%%% FIX %%%%%%%
%d_y=1/3
%d_x=1/5
%%% END FIX %%%

h_y=0.95*d_y;
h_x=h_y*ar;

y=1-d_y:-0.99/n_y:0;
x=0:d_x:1-d_x;

[imY,imX]=meshgrid(y,x);

if textsize==0
    textsize=round(12*4/n_x);
    if textsize<8
        textsize=8;
    end
end
l_t=1/9*textsize/12; % text length
h_t=1/50*textsize/12; % text height
x_t=round((h_x-l_t)/h_x/2*n_j);
x_master=round((h_x-l_t*0.8)/h_x/2*n_j);
y_t=round(h_t*1.3/h_y*n_i);
y_master=round(h_t*6/h_y*n_i);
y_t2=round(h_t*2.5/h_y*n_i);

day=datenum(num2str(ifgname'),'yyyymmdd');
if ~isempty(ifgname2)
    day2=datenum(num2str(ifgname2'),'yyyymmdd');
end
figure
for i=1:n_k
    axes('position',[imX(i),imY(i),h_x,h_y])
    %subplot(6,5,i);
    if heading <90 & heading >-90
        image(flipud(ifgstack(:,:,i)));
    else
        image(fliplr(ifgstack(:,:,i)));
    end
    box on
    axis equal
    axis tight
    set(gca,'yticklabel',[])
    set(gca,'xticklabel',[])
    axis equal
    axis tight
    t=text(x_t,y_t,[datestr(day(i),'dd mmm yyyy')]);
    set(t,'fontweight','bold','color',[1 1 0.996],'fontsize',textsize)
    if exist('day2','var') & length(day2) >= i
        t=text(x_t,y_t2,[datestr(day2(i),'dd mmm yyyy')]);
        set(t,'fontweight','bold','color',[1 1 0.996],'fontsize',textsize)
    end
    if i==master_ix
        t=text(x_master,y_master,'MASTER');
        set(t,'fontweight','bold','color',[1 1 0.996],'fontsize',textsize)
    end        

    %xx=text(0,size(ix,1)*1.15,datestr(day(i)));
    %set(xx,'fontsize',8)

end

colormap(cc)




