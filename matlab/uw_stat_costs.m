function []=uw_stat_costs(unwrap_method,subset_ifg_index);
%UW_STAT_COSTS Find unwrapped solutions using MAP cost functions
%
%   Andy Hooper May 2007
%
%   ======================================================================
%   01/2009 AH: Remove variation due to estimated defo from cost function
%   03/2009 AH: Calculate phase jumps per IFG
%   11/2009 AH: Fix 2D processing
%   ======================================================================


tic
if nargin<1
    unwrap_method='3D';
end

costscale=100;
nshortcycle=200;
maxshort=32000;

fprintf('Unwrapping in space...\n')

uw=load('uw_grid');
ui=load('uw_interp');
ut=load('uw_space_time');

if nargin<2
    subset_ifg_index=[1:size(uw.ph,2)];
end

[nrow,ncol]=size(uw.nzix);

[y,x]=find(uw.nzix);
nzix=find(uw.nzix);
z=[1:uw.n_ps];
 

colix=ui.colix;
rowix=ui.rowix;
Z=ui.Z;


grid_edges=[colix(abs(colix)>0);rowix(abs(rowix)>0)];
n_edges=hist(abs(grid_edges),[1:ui.n_edge])';

if strcmpi(unwrap_method,'2D')
    sigsq_noise=ones(ui.n_edge,1);
    dph_smooth=angle(ut.dph_space);
else
    sigsq_noise=(std(ut.dph_noise,0,2)/2/pi).^2;
    sigsq_defo=(std(ut.dph_space_uw-ut.dph_noise,0,2)/2/pi).^2;
    dph_smooth=ut.dph_space_uw-ut.dph_noise;
end
 
  nostats_ix=find(isnan(sigsq_noise))'; %noise set to nans in uw_sb_unwrap_time
  for i=nostats_ix
    rowix(abs(rowix)==i)=nan;
    colix(abs(colix)==i)=nan;
  end

  sigsq=int16(round(((sigsq_noise)*nshortcycle^2)/costscale.*n_edges)); % weight by number of occurences

  sigsq(sigsq<1)=1; % zero causes snaphu to crash

  rowcost=zeros((nrow-1),ncol*4,'int16');
  colcost=zeros((nrow),(ncol-1)*4,'int16');

  nzrowix=abs(rowix)>0; % 0 = same node, NaN = no stats
  stdgrid=ones(size(rowix),'int16');
  stdgrid(nzrowix)= sigsq(abs(rowix(nzrowix)));
  rowcost(:,2:4:end)= stdgrid; % sigsq

  nzcolix=abs(colix)>0;
  stdgrid=ones(size(colix),'int16');
  stdgrid(nzcolix)= sigsq(abs(colix(nzcolix)));
  colcost(:,2:4:end)= stdgrid; % sigsq

  rowcost(:,3:4:end)= maxshort;% dzmax - value at which cost soars proportional to square of value exceeding dzmax
  colcost(:,3:4:end)= maxshort; % dzmax

  stats_ix=~isnan(rowix); 
  rowcost(:,4:4:end)= int16(stats_ix)*(-1-maxshort)+1; % laycost, -32000 signifies no cost shelf 
  stats_ix=~isnan(colix); 
  colcost(:,4:4:end)= int16(stats_ix)*(-1-maxshort)+1;

ph_uw=zeros(uw.n_ps,uw.n_ifg,'single');
ifguw=zeros(nrow,ncol);
msd=zeros(uw.n_ifg,1);

for i1=subset_ifg_index
    fprintf('\nProcessing IFG %d of %d\n',i1,length(subset_ifg_index));

    offset_cycle=(angle(ut.dph_space(:,i1))-dph_smooth(:,i1))/2/pi;
    offgrid=zeros(size(rowix),'int16');
    offgrid(nzrowix)=round(offset_cycle(abs(rowix(nzrowix))).*sign(rowix(nzrowix))*nshortcycle);
    rowcost(:,1:4:end)= -offgrid; % offset
    offgrid=zeros(size(colix),'int16');
    offgrid(nzcolix)=round(offset_cycle(abs(colix(nzcolix))).*sign(colix(nzcolix))*nshortcycle);
    colcost(:,1:4:end)= offgrid; % offset
    fid=fopen('snaphu.costinfile','w');
    fwrite(fid,rowcost','int16');
    fwrite(fid,colcost','int16');
    fclose(fid);
    ifgw=reshape(uw.ph(Z,i1),nrow,ncol);
    writecpx('snaphu.in',ifgw)
    cmdstr=['!snaphu -d -f $STAMPS/snaphu.conf ',num2str(ncol),' >> snaphu.log'];
    eval(cmdstr)
    fid=fopen('snaphu.out');
    ifguw=fread(fid,[ncol,inf],'float');
    fclose(fid);
    ifguw=ifguw';
    ifg_diff1=ifguw(1:end-1,:)-ifguw(2:end,:);
    ifg_diff1=ifg_diff1(ifg_diff1~=0);
    ifg_diff2=ifguw(:,1:end-1)-ifguw(:,2:end);
    ifg_diff2=ifg_diff2(ifg_diff2~=0);
    msd(i1)=(sum(ifg_diff1.^2)+sum(ifg_diff2.^2))/(length(ifg_diff1)+length(ifg_diff2));
    ph_uw(:,i1)=ifguw(uw.nzix);
end


save('uw_phaseuw','ph_uw','msd')


