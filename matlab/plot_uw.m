function []=plot_uw();
%PLOT_UW plot some figures of intermediate unwrapping steps
%
% Andy Hooper, March 2012

load uw_interp
load uw_grid
load uw_space_time
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
end
set(gcf,'name','Interpolated phase')
subplot(ni,nj,i+1)
imagesc(1,[-pi,pi])
colorbar('westoutside')
cla
axis off

if ~isempty(ifreq_ij)
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
    axis equal
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
    axis equal
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
