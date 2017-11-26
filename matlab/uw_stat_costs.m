function []=uw_stat_costs(unwrap_method,variance,subset_ifg_index);
%UW_STAT_COSTS Find unwrapped solutions using MAP cost functions
%
%   Andy Hooper May 2007
%
%   ======================================================================
%   01/2009 AH: Remove variation due to estimated defo from cost function
%   03/2009 AH: Calculate phase jumps per IFG
%   11/2009 AH: Fix 2D processing
%   01/2012 AH: Create snaphu.conf file internally
%   01/2012 AH: Change noise estimation for 2D method
%   02/2012 AH: Updated for 3D_NEW method
%   03/2013 AH: Add variance option for 2D method
%   08/2014 DB: Suppress external command window output
%   08/2017 AH: Allow for methods that don't set predef_ix
%   ======================================================================


tic
if nargin<1
    unwrap_method='3D';
end

costscale=100;
nshortcycle=200;
maxshort=32000;

fprintf('Unwrapping in space...\n')

uw=load('uw_grid','ph','nzix','pix_size','n_ps','n_ifg');
ui=load('uw_interp');
ut=load('uw_space_time','dph_space_uw','dph_noise','spread','predef_ix');

if nargin<2
    variance=[];
end

if nargin<3
    subset_ifg_index=[1:size(uw.ph,2)];
end

predef_flag='n';
if isfield(ut,'predef_ix') & ~isempty(ut.predef_ix)
    predef_flag='y';
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
    edge_length=sqrt(diff(x(ui.edgs(:,2:3)),[],2).^2+diff(y(ui.edgs(:,2:3)),[],2).^2);
    %sigsq_noise=ones(ui.n_edge,1);
    if uw.pix_size==0
        pix_size=5;  % if we don't know resolution
    else
        pix_size=uw.pix_size;
    end
    if isempty(variance)
        sigsq_noise=zeros(size(edge_length));
    else
        sigsq_noise=variance(ui.edgs(:,2))+variance(ui.edgs(:,3));
    end
    sigsq_aps=(2*pi)^2; % fixed for now as one fringe
    aps_range=20000; % fixed for now as 20 km
    sigsq_noise=sigsq_noise+sigsq_aps*(1-exp(-edge_length*pix_size*3/aps_range)); % cov of dph=C11+C22-2*C12 (assume APS only contributor)
    sigsq_noise=sigsq_noise/10; % scale it to keep in reasonable range
    dph_smooth=ut.dph_space_uw;
else
    sigsq_noise=(std(ut.dph_noise,0,2)/2/pi).^2;
    %sigsq_defo=(std(ut.dph_space_uw-ut.dph_noise,0,2)/2/pi).^2;
    dph_smooth=ut.dph_space_uw-ut.dph_noise;
end
ut=rmfield(ut,'dph_noise');

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
  rowstdgrid=ones(size(rowix),'int16');

  nzcolix=abs(colix)>0;
  colstdgrid=ones(size(colix),'int16');
  
  rowcost(:,3:4:end)= maxshort;% dzmax - value at which cost soars proportional to square of value exceeding dzmax
  colcost(:,3:4:end)= maxshort; % dzmax

  stats_ix=~isnan(rowix); 
  rowcost(:,4:4:end)= int16(stats_ix)*(-1-maxshort)+1; % laycost, -32000 signifies no cost shelf 
  stats_ix=~isnan(colix); 
  colcost(:,4:4:end)= int16(stats_ix)*(-1-maxshort)+1;

ph_uw=zeros(uw.n_ps,uw.n_ifg,'single');
ifguw=zeros(nrow,ncol);
msd=zeros(uw.n_ifg,1);

fid=fopen('snaphu.conf','w');
fprintf(fid,'INFILE  snaphu.in\n');
fprintf(fid,'OUTFILE snaphu.out\n');
fprintf(fid,'COSTINFILE snaphu.costinfile\n');
fprintf(fid,'STATCOSTMODE  DEFO\n');
fprintf(fid,'INFILEFORMAT  COMPLEX_DATA\n');
fprintf(fid,'OUTFILEFORMAT FLOAT_DATA\n');
fclose(fid);


for i1=subset_ifg_index
    fprintf('   Processing IFG %d of %d\n',i1,length(subset_ifg_index));
    spread=full(ut.spread(:,i1));
    spread=int16(round((abs(spread)*nshortcycle^2)/6/costscale.*repmat(n_edges,1,size(spread,2))));
    sigsqtot=sigsq+spread;
    if predef_flag=='y'
        sigsqtot(ut.predef_ix(:,i1))=1;
    end
    rowstdgrid(nzrowix)= sigsqtot(abs(rowix(nzrowix)));
    rowcost(:,2:4:end)= rowstdgrid; % sigsq
    colstdgrid(nzcolix)= sigsqtot(abs(colix(nzcolix)));
    colcost(:,2:4:end)= colstdgrid; % sigsq
    
%    offset_cycle=(angle(ut.dph_space(:,i1))-dph_smooth(:,i1))/2/pi;
    offset_cycle=(angle(exp(1i*ut.dph_space_uw(:,i1)))-dph_smooth(:,i1))/2/pi;
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
    cmdstr=['snaphu -d -f snaphu.conf ',num2str(ncol),' >& snaphu.log'];
    
    
    [a,b] =system(cmdstr);
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


