function []=ps_merge_patches(psver)
%PS_MERGE_PATCHES merge patches
%
%   Andy Hooper, September 2006
%  
%   ======================================================================
%   01/2007 AH: scla and scn added
%   01/2007 AH: patch.list added
%   06/2009 AH: file existence checks to search current directory only
%   06/2009 AH: move mean amplitude merge to end to save memory
%   06/2009 AH: only save scla and scn when present to save memory
%   09/2009 AH: add option to resample to coarser sampling 
%   09/2009 AH: reduce memory needs further
%   11/2009 AH: ensure mean amplitude width is always correct
%   06/2010 AH: estimate weights directly from residuals
%   10/2010 DB: Fix when patch_noover does not have PS (resampling)
%   10/2010 DB: Fix when PS all rejected when sum weight<min_weight (resampling)
%   06/2010 AH: Move mean amplitude merging to ps_load_mean_amp.m 
%   02/2011 DB: Fix dimension for min computation of ps.xy 
%   09/2015 AH: Delete previous merged amplitude files
%   06/2017 DB: Include stamps save for large variables
%   10/2017 DB: If inc angle file is present also merge it.
%   ======================================================================

logit;
fprintf('Merging patches...\n')


if nargin < 1
  psver=2;
end

small_baseline_flag=getparm('small_baseline_flag');
grid_size=getparm('merge_resample_size',1);
merge_stdev=getparm('merge_standard_dev',1);
phase_accuracy=10*pi/180; % gives minimum possible accuracy for a pixel 
min_weight=1/merge_stdev^2;
randn('state',1001);
max_coh=abs(sum(exp(j*randn(1000,1)*phase_accuracy)))/1000;

psname=['ps',num2str(psver)];
phname=['ph',num2str(psver)];
rcname=['rc',num2str(psver)];
pmname=['pm',num2str(psver)];
phuwname=['phuw',num2str(psver)];
sclaname=['scla',num2str(psver)];
sclasbname=['scla_sb',num2str(psver)];
scnname=['scn',num2str(psver)];
bpname=['bp',num2str(psver)];
laname=['la',num2str(psver)];
incname=['inc',num2str(psver)];
hgtname=['hgt',num2str(psver)];

if exist('./patch.list','file')
    dirname=struct;
    fid=fopen('patch.list','r');
    i=0;
    while feof(fid)==0
        i=i+1;
        dirname(i).name=fgetl(fid);
    end
    fclose(fid);
else
    dirname=dir('PATCH_*');
end

n_patch=length(dirname);
remove_ix=logical([]);
ij=zeros(0,2);
lonlat=zeros(0,2);
ph=zeros(0,0);
ph_rc=zeros(0,0);
ph_reref=zeros(0,0);
ph_uw=zeros(0,0);
ph_patch=zeros(0,0);
ph_res=zeros(0,0);
ph_scla=zeros(0,0,'single');
ph_scla_sb=zeros(0,0,'single');
ph_scn_master=zeros(0,0);
ph_scn_slave=zeros(0,0);
K_ps=zeros(0,0);
C_ps=zeros(0,0);
coh_ps=zeros(0,0);
K_ps_uw=zeros(0,0,'single');
K_ps_uw_sb=zeros(0,0,'single');
C_ps_uw=zeros(0,0,'single');
C_ps_uw_sb=zeros(0,0,'single');
bperp_mat=zeros(0,0,'single');
la=zeros(0,0);
inc=zeros(0,0);
hgt=zeros(0,0);
amp=zeros(0,0,'single');

for i=1:n_patch
  if ~isempty(dirname(i).name)
    fprintf('   Processing %s\n',dirname(i).name)
    cd(dirname(i).name);
    ps=load(psname);
    n_ifg=ps.n_ifg;
    if isfield('ps','n_image')
        n_image=ps.n_image;
    else
        n_image=ps.n_ifg;
    end
    
    patch.ij=load('patch_noover.in');
    ix=(ps.ij(:,2)>=patch.ij(3)-1 & ps.ij(:,2)<=patch.ij(4)-1 & ps.ij(:,3)>=patch.ij(1)-1 & ps.ij(:,3)<=patch.ij(2)-1);
    if sum(ix)==0
   	ix_no_ps =1;	% no PS left afer removing overlapping patches
    else
	ix_no_ps=0;
    end 

    if grid_size==0
      [C,IA,IB]=intersect(ps.ij(ix,2:3),ij,'rows'); 
      remove_ix=[remove_ix;IB]; % now more reliable values for these pixels
      [C,IA,IB]=intersect(ps.ij(:,2:3),ij,'rows');
      ix_ex=true(ps.n_ps,1);
      ix_ex(IA)=0; % exclusive index (non-intersecting)
      ix(ix_ex)=1; % keep those in patch proper + those outside not already kept
    elseif grid_size ~=0 && ix_no_ps~=1
      clear g_ij
      xy_min=min(ps.xy(ix,:),1);
      g_ij(:,1)=ceil((ps.xy(ix,3)-xy_min(3)+1e-9)/grid_size);
      g_ij(:,2)=ceil((ps.xy(ix,2)-xy_min(2)+1e-9)/grid_size);
      n_i=max(g_ij(:,1));
      n_j=max(g_ij(:,2));
      [g_ij,I,g_ix]=unique(g_ij,'rows');
      [g_ix,sort_ix]=sort(g_ix);
      ix=find(ix);
      ix=ix(sort_ix);
      pm=load(pmname,'ph_res','coh_ps','C_ps');
      pm.ph_res=angle(exp(j*(pm.ph_res-repmat(pm.C_ps,1,size(pm.ph_res,2))))); % centralise about zero
      if small_baseline_flag~='y'
          pm.ph_res=[pm.ph_res,pm.C_ps]; %  include master noise too
      end
      sigsq_noise=var([pm.ph_res],0,2); 
      coh_ps_all=abs(sum(exp(j*[pm.ph_res]),2))/n_ifg;
      coh_ps_all(coh_ps_all>max_coh)=max_coh; % % prevent unrealistic weights
      sigsq_noise(sigsq_noise<phase_accuracy^2)=phase_accuracy^2; % prevent unrealistic weights
      ps_weight=1./sigsq_noise(ix);
      ps_snr=1./(1./coh_ps_all(ix).^2-1);
      clear pm

      l_ix=[find(diff(g_ix));size(g_ix,1)];
      f_ix=[1;l_ix(1:end-1)+1];
      n_ps_g=size(f_ix,1); 

      weightsave=zeros(n_ps_g,1);
      for i=1:n_ps_g
          weights=ps_weight(f_ix(i):l_ix(i));
          weightsum=sum(weights);
          weightsave(i)=weightsum;
          if weightsave(i)<min_weight
              ix(f_ix(i):l_ix(i))=0;
          end
      end
      g_ix=g_ix(ix>0);
      
      if isempty(g_ix)==1
	  ix_no_ps=1;				% Remaining PS are rejected because to sum of weights is smaller than threshold min_weight		
      end

      l_ix=[find(diff(g_ix));size(g_ix,1)];
      f_ix=[1;l_ix(1:end-1)+1];
      ps_weight=ps_weight(ix>0);
      ps_snr=ps_snr(ix>0);
      ix=ix(ix>0);
      n_ps_g=size(f_ix,1); 
      n_ps=length(ix);
    end
    
    if grid_size==0
      ij=[ij;ps.ij(ix,2:3)];
      lonlat=[lonlat;ps.lonlat(ix,:)];
    elseif grid_size ~=0 && ix_no_ps~=1
      ij_g=zeros(n_ps_g,2);
      lonlat_g=zeros(n_ps_g,2);
      ps.ij=ps.ij(ix,:);
      ps.lonlat=ps.lonlat(ix,:);
      for i=1:n_ps_g
        weights=repmat(ps_weight(f_ix(i):l_ix(i)),1,2);
        ij_g(i,:)=round(sum(ps.ij(f_ix(i):l_ix(i),2:3).*weights,1)/sum(weights(:,1)));
        lonlat_g(i,:)=sum(ps.lonlat(f_ix(i):l_ix(i),:).*weights,1)/sum(weights(:,1));
      end
      ij=[ij;ij_g];
      lonlat=[lonlat;lonlat_g];
    end


    if exist(['./',phname,'.mat'],'file')
        phin=load(phname);
        ph_w=phin.ph;
        clear phin;
    elseif isfield(ps,'ph')
        ph_w=ps.ph;
    end
    if exist('ph_w','var')
      if grid_size==0
        ph=[ph;ph_w(ix,:)];
      elseif grid_size ~=0 && ix_no_ps~=1
        ph_w=ph_w(ix,:);
        ph_g=zeros(n_ps_g,n_ifg);
        for i=1:n_ps_g
          weights=repmat(ps_snr(f_ix(i):l_ix(i)),1,n_ifg);
          ph_g(i,:)=sum(ph_w(f_ix(i):l_ix(i),:).*weights,1);
        end
        ph=[ph;ph_g];
        clear ph_g
      end
      clear ph_w
    end

    rc=load(rcname);
    if grid_size==0
      ph_rc=[ph_rc;rc.ph_rc(ix,:)];
      if ~strcmpi(small_baseline_flag,'y')
        ph_reref=[ph_reref;rc.ph_reref(ix,:)];
      end
    elseif grid_size ~=0 && ix_no_ps~=1
      rc.ph_rc=rc.ph_rc(ix,:);
      ph_g=zeros(n_ps_g,n_ifg);
      if ~strcmpi(small_baseline_flag,'y')
        rc.ph_reref=rc.ph_reref(ix,:);
        ph_reref_g=zeros(n_ps_g,n_ifg);
      end
      for i=1:n_ps_g
        weights=repmat(ps_snr(f_ix(i):l_ix(i)),1,n_ifg);
        ph_g(i,:)=sum(rc.ph_rc(f_ix(i):l_ix(i),:).*weights,1);
        if ~strcmpi(small_baseline_flag,'y')
          ph_reref_g(i,:)=sum(rc.ph_reref(f_ix(i):l_ix(i),:).*weights,1);
        end
      end
      ph_rc=[ph_rc;ph_g];
      clear ph_g
      if ~strcmpi(small_baseline_flag,'y')
        ph_reref=[ph_reref;ph_reref_g];
        clear ph_reref_g
      end
    end
    clear rc 
    
    pm=load(pmname);
    if grid_size==0
      ph_patch=[ph_patch;pm.ph_patch(ix,:)];
      if isfield(pm,'ph_res')
        ph_res=[ph_res;pm.ph_res(ix,:)];
      end
      if isfield(pm,'K_ps')
        K_ps=[K_ps;pm.K_ps(ix,:)];
      end
      if isfield(pm,'C_ps')
        C_ps=[C_ps;pm.C_ps(ix,:)];
      end
      if isfield(pm,'coh_ps')
        coh_ps=[coh_ps;pm.coh_ps(ix,:)];
      end
    elseif grid_size ~=0 && ix_no_ps~=1
      pm.ph_patch=pm.ph_patch(ix,:);
      ph_g=zeros(n_ps_g,size(pm.ph_patch,2));
      if isfield(pm,'ph_res')
        pm.ph_res=pm.ph_res(ix,:);
        ph_res_g=ph_g;
      end
      if isfield(pm,'K_ps')
        pm.K_ps=pm.K_ps(ix,:);
        K_ps_g=zeros(n_ps_g,1);
      end
      if isfield(pm,'C_ps')
        pm.C_ps=pm.C_ps(ix,:);
        C_ps_g=zeros(n_ps_g,1);
      end
      if isfield(pm,'coh_ps')
        pm.coh_ps=pm.coh_ps(ix,:);
        coh_ps_g=zeros(n_ps_g,1);
      end
      for i=1:n_ps_g
        weights=repmat(ps_snr(f_ix(i):l_ix(i)),1,size(ph_g,2));
        ph_g(i,:)=sum(pm.ph_patch(f_ix(i):l_ix(i),:).*weights,1);
        if isfield(pm,'ph_res')
          ph_res_g(i,:)=sum(pm.ph_res(f_ix(i):l_ix(i),:).*weights,1);
        end
        if isfield(pm,'coh_ps')
          snr=sqrt(sum(weights(:,1).^2,1));
          coh_ps_g(i)=sqrt(1./(1+1./snr));
        end
        weights=ps_weight(f_ix(i):l_ix(i));
        if isfield(pm,'K_ps')
          K_ps_g(i)=sum(pm.K_ps(f_ix(i):l_ix(i),:).*weights,1)./sum(weights,1);
        end
        if isfield(pm,'C_ps')
          C_ps_g(i)=sum(pm.C_ps(f_ix(i):l_ix(i),:).*weights,1)./sum(weights,1);
        end

	if sum(sum(isnan(C_ps_g)))>0 || sum(sum(isnan(weights)))>0 || sum(sum(isnan(ph_g)))>0 ||  sum(sum(isnan(K_ps_g)))>0 ||  sum(sum(isnan(coh_ps_g)))>0 ||  sum(sum(isnan(snr)))>0
keyboard
end 

      end
      ph_patch=[ph_patch;ph_g];
      clear ph_g
      if isfield(pm,'ph_res')
        ph_res=[ph_res;ph_res_g];
        clear ph_res_g
      end
      if isfield(pm,'K_ps')
        K_ps=[K_ps;K_ps_g];
        clear K_ps_g
      end
      if isfield(pm,'C_ps')
        C_ps=[C_ps;C_ps_g];
        clear C_ps_g
      end
      if isfield(pm,'coh_ps')
        coh_ps=[coh_ps;coh_ps_g];
        clear coh_ps_g
      end
    end
    clear pm

    bp=load(bpname);
    if grid_size==0
      bperp_mat=[bperp_mat;bp.bperp_mat(ix,:)];
    elseif grid_size ~=0 && ix_no_ps~=1
      bperp_g=zeros(n_ps_g,size(bp.bperp_mat,2));
      bp.bperp_mat=bp.bperp_mat(ix,:);
      for i=1:n_ps_g
        weights=repmat(ps_weight(f_ix(i):l_ix(i)),1,size(bperp_g,2));
        weights(weights==0)=1e-9; % pixels with zero phase cause this problem
        bperp_g(i,:)=sum(bp.bperp_mat(f_ix(i):l_ix(i),:).*weights,1)/sum(weights(:,1));
      end
      bperp_mat=[bperp_mat;bperp_g];
      clear bperp_g
    end
    clear bp

    if exist(['./',laname,'.mat'],'file')
      lain=load(laname);
      if grid_size==0
        la=[la;lain.la(ix,:)];
      elseif grid_size ~=0 && ix_no_ps~=1
        la_g=zeros(n_ps_g,1);
        lain.la=lain.la(ix,:);
        for i=1:n_ps_g
          weights=ps_weight(f_ix(i):l_ix(i));
          la_g(i)=sum(lain.la(f_ix(i):l_ix(i)).*weights,1)/sum(weights(:,1));
        end  
        la=[la;la_g];
        clear la_g
      end
      clear lain
    end

     if exist(['./',incname,'.mat'],'file')
      incin=load(incname);
      if grid_size==0
        inc=[inc;incin.inc(ix,:)];
      elseif grid_size ~=0 && ix_no_ps~=1
        inc_g=zeros(n_ps_g,1);
        incin.inc=incin.inc(ix,:);
        for i=1:n_ps_g
          weights=ps_weight(f_ix(i):l_ix(i));
          inc_g(i)=sum(incin.inc(f_ix(i):l_ix(i)).*weights,1)/sum(weights(:,1));
        end  
        inc=[inc;inc_g];
        clear inc_g
      end
      clear incin
    end
    
    
    
    if exist(['./',hgtname,'.mat'],'file')
      hgtin=load(hgtname);
      if grid_size==0
        hgt=[hgt;hgtin.hgt(ix,:)];
      elseif grid_size ~=0 && ix_no_ps~=1
        hgt_g=zeros(n_ps_g,1);
        hgtin.hgt=hgtin.hgt(ix,:);
        for i=1:n_ps_g
          weights=ps_weight(f_ix(i):l_ix(i));
          hgt_g(i)=sum(hgtin.hgt(f_ix(i):l_ix(i)).*weights,1)/sum(weights(:,1));
        end  
        hgt=[hgt;hgt_g];
        clear hgt_g
      end
      clear hgtin
    end

    if grid_size==0
    %overlap_ij=ps.ij(~ix,:);
      if exist(['./',phuwname,'.mat'],'file')
        phuw=load(phuwname);
        if ~isempty(C)
            ph_uw_diff=mean(phuw.ph_uw(IA,:)-ph_uw(IB,:),1);
            if ~strcmpi(small_baseline_flag,'y')
                ph_uw_diff=round(ph_uw_diff/2/pi)*2*pi; % round to nearest 2 pi
            end
        else
            ph_uw_diff=zeros(1,size(phuw.ph_uw,2));
        end
        ph_uw=[ph_uw;phuw.ph_uw(ix,:)-repmat(ph_uw_diff,sum(ix),1)];
        clear phuw
      else
        ph_uw=[ph_uw;zeros(sum(ix),n_image,'single')];
      end
    
      if exist(['./',sclaname,'.mat'],'file')
        scla=load(sclaname);
        if ~isempty(C)
            ph_scla_diff=mean(scla.ph_scla(IA,:)-ph_scla(IB,:));
            K_ps_diff=mean(scla.K_ps_uw(IA,:)-K_ps_uw(IB,:));
            C_ps_diff=mean(scla.C_ps_uw(IA,:)-C_ps_uw(IB,:));
        else
            ph_scla_diff=zeros(1,size(scla.ph_scla,2));
            K_ps_diff=0;
            C_ps_diff=0;
        end
        ph_scla=[ph_scla;scla.ph_scla(ix,:)-repmat(ph_scla_diff,sum(ix),1)];
        K_ps_uw=[K_ps_uw;scla.K_ps_uw(ix,:)-repmat(K_ps_diff,sum(ix),1)];
        C_ps_uw=[C_ps_uw;scla.C_ps_uw(ix,:)-repmat(C_ps_diff,sum(ix),1)];
        clear scla
      end

      if strcmpi(small_baseline_flag,'y')
       if exist(['./',sclasbname,'.mat'],'file')
        sclasb=load(sclasbname);
        ph_scla_diff=mean(sclasb.ph_scla(IA,:)-ph_scla_sb(IB,:));
        K_ps_diff=mean(sclasb.K_ps_uw(IA,:)-K_ps_uw_sb(IB,:));
        C_ps_diff=mean(sclasb.C_ps_uw(IA,:)-C_ps_uw_sb(IB,:));
        ph_scla_sb=[ph_scla_sb;sclasb.ph_scla(ix,:)-repmat(ph_scla_diff,sum(ix),1)];
        K_ps_uw_sb=[K_ps_uw_sb;sclasb.K_ps_uw(ix,:)-repmat(K_ps_diff,sum(ix),1)];
        C_ps_uw_sb=[C_ps_uw_sb;sclasb.C_ps_uw(ix,:)-repmat(C_ps_diff,sum(ix),1)];
         clear sclasb
       end
      end

    
    
      if exist(['./',scnname,'.mat'],'file')
        scn=load(scnname);

        if ~isempty(C)
            ph_scn_diff=mean(scn.ph_scn_slave(IA,:)-ph_scn_slave(IB,:));
        else
            ph_scn_diff=zeros(1,size(scn.ph_scn_slave,2));
        end
        ph_scn_slave=[ph_scn_slave;scn.ph_scn_slave(ix,:)-repmat(ph_scn_diff,sum(ix),1)];
        clear scn
      end
    end
    cd ..
  end
end

ps_new=ps;
n_ps_orig=size(ij,1); % before duplicates removed
keep_ix=true(n_ps_orig,1);
keep_ix(remove_ix)=0;
lonlat_save=lonlat;
coh_ps_weed=coh_ps(keep_ix);
lonlat=lonlat(keep_ix,:);

%%%%%%%%% Some non-adjacent pixels are allocated the same lon/lat by DORIS. %%%%%%%%%
%%%%%%%%% If duplicates occur, the pixel with the highest coherence is kept.%%%%%%%%% 
[dummy,I]=unique(lonlat,'rows');
dups=setxor(I,[1:size(lonlat,1)]'); % pixels with duplicate lon/lat
keep_ix_num=find(keep_ix);

for i=1:length(dups)
    dups_ix_weed=find(lonlat(:,1)==lonlat(dups(i),1)&lonlat(:,2)==lonlat(dups(i),2));
    dups_ix=keep_ix_num(dups_ix_weed);
    [dummy,I]=max(coh_ps_weed(dups_ix_weed));
    keep_ix(dups_ix([1:end]~=I))=0; % drop dups with lowest coh
end

if ~isempty(dups)
    lonlat=lonlat_save(keep_ix,:);
    fprintf('   %d pixel with duplicate lon/lat dropped\n\n',length(dups)')
end
clear lonlat_save
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

ll0=(max(lonlat)+min(lonlat))/2;
xy=llh2local(lonlat',ll0)*1000;
xy=xy';
sort_x=sortrows(xy,1);
sort_y=sortrows(xy,2);
n_pc=round(size(xy,1)*0.001);
bl=mean(sort_x(1:n_pc,:)); % bottom left corner
tr=mean(sort_x(end-n_pc:end,:)); % top right corner
br=mean(sort_y(1:n_pc,:)); % bottom right  corner
tl=mean(sort_y(end-n_pc:end,:)); % top left corner

heading=getparm('heading');
if isempty(heading);
    heading=0;
end
theta=(180-heading)*pi/180;
if theta>pi
    theta=theta-2*pi;
end

rotm=[cos(theta),sin(theta); -sin(theta),cos(theta)];
xy=xy';
xynew=rotm*xy; % rotate so that scene axes approx align with x=0 and y=0
if max(xynew(1,:))-min(xynew(1,:))<max(xy(1,:))-min(xy(1,:)) &...
   max(xynew(2,:))-min(xynew(2,:))<max(xy(2,:))-min(xy(2,:))
    xy=xynew; % check that rotation is an improvement
    disp(['   Rotating xy by ',num2str(theta*180/pi),' degrees']);
end
clear xynew

xy=single(xy');
%xy(weightsave<min_weight,:)=1e9;
[xy_sort,sort_ix]=sortrows(xy,[2,1]); % sort in ascending y order
%sort_ix=sort_ix(xy_sort(:,1)<1e9);



xy=xy(sort_ix,:);
xy=[[1:size(xy,1)]',xy];
xy(:,2:3)=round(xy(:,2:3)*1000)/1000; % round to mm
lonlat=lonlat(sort_ix,:);

all_ix=[1:size(ij,1)]';
keep_ix=all_ix(keep_ix);
sort_ix=keep_ix(sort_ix);

n_ps=length(sort_ix);
fprintf('   Writing merged dataset (contains %d pixels)\n',n_ps)

ij=ij(sort_ix,:);

ph_rc=ph_rc(sort_ix,:);
ph_rc(ph_rc~=0)=ph_rc(ph_rc~=0)./abs(ph_rc(ph_rc~=0));
if ~strcmpi(small_baseline_flag,'y')
    ph_reref=ph_reref(sort_ix,:);
end
stamps_save(rcname,ph_rc,ph_reref);
clear ph_rc ph_reref

if size(ph_uw,1)==n_ps_orig
    ph_uw=ph_uw(sort_ix,:);
    stamps_save(phuwname,ph_uw);
end
clear ph_uw

ph_patch=ph_patch(sort_ix,:);
if size(ph_res,1)==n_ps_orig
    ph_res=ph_res(sort_ix,:);
else
    ph_res=[];
end
if size(K_ps,1)==n_ps_orig
    K_ps=K_ps(sort_ix,:);
else
    K_ps=[];
end
if size(C_ps,1)==n_ps_orig
    C_ps=C_ps(sort_ix,:);
else
    C_ps=[];
end
if size(coh_ps,1)==n_ps_orig
    coh_ps=coh_ps(sort_ix,:);
else
    coh_ps=[];
end
stamps_save(pmname,ph_patch,ph_res,K_ps,C_ps,coh_ps);
clear ph_patch ph_res K_ps C_ps coh_ps

if size(ph_scla,1)==n_ps
    ph_scla=ph_scla(sort_ix,:);
    K_ps_uw=K_ps_uw(sort_ix,:);
    C_ps_uw=C_ps_uw(sort_ix,:);
    stamps_save(sclaname,ph_scla,K_ps_uw,C_ps_uw);
end
clear ph_scla K_ps_uw C_ps_uw

if strcmpi(small_baseline_flag,'y') & size(ph_scla_sb,1)==n_ps
  ph_scla=ph_scla_sb(sort_ix,:);
  K_ps_uw=K_ps_uw_sb(sort_ix,:);
  C_ps_uw=C_ps_uw_sb(sort_ix,:);
  stamps_save(sclasbname,ph_scla,K_ps_uw,C_ps_uw);
  clear ph_scla K_ps_uw C_ps_uw 
end
clear ph_scla_sb K_ps_uw_sb C_ps_uw_sb

if size(ph_scn_slave,1)==n_ps
  %ph_scn_master=ph_scn_master(sort_ix,:);
  ph_scn_slave=ph_scn_slave(sort_ix,:);
  stamps_save(scnname,ph_scn_slave);
end
clear ph_scn_slave

if size(ph,1)==n_ps_orig
    ph=ph(sort_ix,:);
else
    ph=[];
end
stamps_save(phname,ph);
clear ph


if size(la,1)==n_ps_orig
    la=la(sort_ix,:);
else
    la=[];
end
stamps_save(laname,la);
clear la

if size(inc,1)==n_ps_orig
    inc=inc(sort_ix,:);
else
    inc=[];
end
stamps_save(incname,inc);
clear inc



if size(hgt,1)==n_ps_orig
    hgt=hgt(sort_ix,:);
else
    hgt=[];
end
stamps_save(hgtname,hgt);
clear hgt

bperp_mat=bperp_mat(sort_ix,:);
stamps_save(bpname,bperp_mat);
clear bperp_mat

ps_new.n_ps=n_ps;
ps_new.ij=[[1:n_ps]',ij];
ps_new.xy=xy;
ps_new.lonlat=lonlat;
save(psname,'-struct','ps_new');
clear ps_new

save psver psver

vars=who;
vars=setxor(vars,{'n_patch';'dirname'});
clear(vars{:})

if exist('./mean_amp.flt','file')  
    delete('./mean_amp.flt') 
end

if exist('./amp_mean.mat','file')  
    delete('./amp_mean.mat') 
end


logit(1);
