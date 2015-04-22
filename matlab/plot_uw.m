function []=plot_uw(ifg_list);
%PLOT_UW plot some figures of intermediate unwrapping steps
%   PLOT_UW(IFG_LIST) - default is to plot the first
%
% Andy Hooper, March 2012
% modifications:
% DB    11/2014         Show the interferogram dates when plotted  

if nargin<1
    ifg_list=1;
end

load uw_interp
load uw_grid
load ps2

drop_ifg = getparm('drop_ifg');
ifg_keep = [1:n_ifg]';
ifg_keep(drop_ifg)=[];
ifg_keep = ifg_keep(ifg_list);
if strcmpi(getparm('small_baseline_flag'),'y')
    ifg_str = [repmat('ifg ',length(ifg_keep),1)  num2str(ifg_keep) repmat(': ',length(ifg_keep),1) datestr(ifgday(ifg_keep,1)) repmat(' till ',length(ifg_keep),1)  datestr(ifgday(ifg_keep,2))];
else
    ifg_str = '';
end
ut=load('uw_space_time')

ph=ph(:,ifg_list);

n_ifg=size(ph,2);
ni=floor(sqrt(n_ifg+1));
nj=ceil(sqrt(n_ifg+1));
if ni*nj<n_ifg+1
    ni=ni+1;
end
[nrow,ncol]=size(Z);

figure
for i=1:n_ifg
     subplot(ni,nj,i)
     ifgw=ones(nrow,ncol)*(pi+pi/31.9);
     ifgw(nzix)=angle(ph(:,i));
     imagesc(ifgw,[-pi,pi+pi/31.9])
     axis off
     axis equal
     axis xy
     if ~isempty(ifg_str)
        title(ifg_str(i,:),'fontsize',15)
     end
end
colormap([jet(64);0,0,0]);
subplot(ni,nj,i+1)
imagesc(1,[-pi,pi])
colorbar('westoutside')
cla
axis off
set(gcf,'name','Resampled (Filtered) Phase')
 
figure
for i=1:n_ifg
    subplot(ni,nj,i)
    ifgw=angle(reshape(ph(Z,i),nrow,ncol));
    imagesc(ifgw,[-pi,pi])
    axis off
    axis equal
    axis xy
    if ~isempty(ifg_str)
        title(ifg_str(i,:),'fontsize',15)
     end
end
set(gcf,'name','Interpolated phase')
subplot(ni,nj,i+1)
imagesc(1,[-pi,pi])
colorbar('westoutside')
cla
axis off

load uw_space_time ifreq_ij jfreq_ij

if ~isempty(ifreq_ij)
    ifreq_ij=ifreq_ij(:,n_ifg);
    jfreq_ij=jfreq_ij(:,n_ifg);
    figure
    for i=1:n_ifg
        subplot(ni,nj,i)
        ifgw=zeros(nrow,ncol);
        clim=[min(ifreq_ij(:)),max(ifreq_ij(:))];
        clim(2)=clim(2)+diff(clim)/63.9;
        ifgw(nzix)=(ifreq_ij(:,i));
        ifgw(~nzix)=clim(2);
        imagesc(ifgw,clim);
        axis off
    %    axis equal
        axis xy
    end
    colormap([jet(64);0,0,0]);
    set(gcf,'name','North Phase Gradient')
    subplot(ni,nj,i+1)
    imagesc(1,clim)
    colorbar('westoutside')
    cla
    axis off
end

if ~isempty(jfreq_ij)
figure
for i=1:n_ifg
    subplot(ni,nj,i)
    ifgw=zeros(nrow,ncol);
    clim=[min(jfreq_ij(:)),max(jfreq_ij(:))];
    clim(2)=clim(2)+diff(clim)/63.9;
    ifgw(nzix)=(jfreq_ij(:,i));
    ifgw(~nzix)=clim(2);
    imagesc(ifgw,clim);
    axis off
%    axis equal
    axis xy
end
colormap([jet(64);0,0,0]);
set(gcf,'name','East Phase Gradient')
subplot(ni,nj,i+1)
imagesc(1,clim)
colorbar('westoutside')
cla
axis off
end

figure
ifgw=ones(nrow,ncol)*(pi+pi/31.9);
ifgw(nzix)=angle(ph(:,1));
imagesc([min(xy(:,2)),max(xy(:,2))],[min(xy(:,3)),max(xy(:,3))],ifgw,[-pi,pi+pi/31.9])
colormap([jet(64);0,0,0]);

plot_edges(edges,xy(:,2),xy(:,3),'w');
if isfield(ut,'shaky_ix')
    plot_edges(edges(ut.shaky_ix,:),xy(:,2),xy(:,3),'r');
end

axis off
axis xy


set(gcf,'name','Arcs')

end

