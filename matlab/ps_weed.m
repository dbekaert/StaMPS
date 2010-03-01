function []=ps_weed(all_da_flag,no_weed_adjacent,no_weed_noisy)
%PS_WEED weeds out neighboring PS and save those kept to new version
%   PS_weed(all_da_flag,no_weed_adjacent,no_weed_noisy)
%
%   Andy Hooper, June 2006
%
%   =================================================================
%   09/2006 AH: create all workspace files directly
%   09/2006 AH: drop noisy pixels
%   09/2006 AH: add small baselines 
%   01/2007 AH: drop pixels with duplicate lon/lat 
%   05/2007 AH: optional weeding of pixels with zero elevation added 
%   03/2009 AH: delete scla_smooth mat files
%   02/2010 AH: change smoothing of arcs to time domain
%   02/2010 AH: option to threshold on max arc noise in any ifg
%   =================================================================

fprintf('Weeding selected pixels...\n')

if nargin<1
    all_da_flag=0;
end

%weed_alpha=getparm('weed_alpha',1);
time_win=getparm('weed_time_win',1);
weed_standard_dev=getparm('weed_standard_dev',1);
weed_max_noise=getparm('weed_max_noise',1);
weed_zero_elevation=getparm('weed_zero_elevation',1);
small_baseline_flag=getparm('small_baseline_flag',1);

if nargin<2
    no_weed_adjacent=0;
end
if nargin<3
    if weed_standard_dev>=pi && weed_max_noise>=pi
        no_weed_noisy=1;
    else
        no_weed_noisy=0;
    end
end

load psver

psname=['ps',num2str(psver)];
pmname=['pm',num2str(psver)];
phname=['ph',num2str(psver)];
selectname=['select',num2str(psver)];
hgtname=['hgt',num2str(psver),'.mat'];
laname=['la',num2str(psver),'.mat'];
bpname=['bp',num2str(psver),'.mat'];

psothername=['ps_other'];
psothername=['pm_other'];
selectothername=['select_other'];
hgtothername=['hgt_other'];
laothername=['la_other'];
bpothername=['bp_other'];

ps=load(psname);
sl=load(selectname);

if exist([phname,'.mat'],'file')
    phin=load(phname);
    ph=phin.ph;
    clear phin
else
    ph=ps.ph;
end

day=ps.day;
bperp=ps.bperp;
master_day=ps.master_day;

if isfield(sl,'keep_ix')
    ix2=sl.ix(sl.keep_ix);
    K_ps2=sl.K_ps2(sl.keep_ix);
    C_ps2=sl.C_ps2(sl.keep_ix);
    coh_ps2=sl.coh_ps2(sl.keep_ix);
else
    ix2=sl.ix2;
    K_ps2=sl.K_ps2;
    C_ps2=sl.C_ps2;
    coh_ps2=sl.coh_ps2;
end

ij2=ps.ij(ix2,:);
xy2=ps.xy(ix2,:);
ph2=ph(ix2,:);
lonlat2=ps.lonlat(ix2,:);

pm=load(pmname);
ph_patch2=pm.ph_patch(ix2,:); % use original patch phase, with PS left in
if isfield(sl,'ph_res2')
    ph_res2=sl.ph_res2(sl.keep_ix,:);
else
    ph_res2=[];
end
clear pm

clear sl

clear ph
if isfield(ps,'ph')
    ps=rmfield(ps,'ph');
end
ps=rmfield(ps,{'xy','ij','lonlat','sort_ix'});

if all_da_flag~=0
    pso=load(psothername);
    slo=load(selectothername);
    ix_other=slo.ix_other;
    n_ps_other=sum(ix_other);
    K_ps_other2=pso.K_ps_other(ix_other);
    C_ps_other2=pso.C_ps_other(ix_other);
    coh_ps_other2=pso.coh_ps_other(ix_other);
    ph_res_other2=pso.ph_res_other(ix_other,:);
    ij2=[ij2;pso.ij_other(ix_other,:)];
    xy2=[xy2;pso.xy_other(ix_other,:)];
    ph2=[ph2;pso.ph_other(ix_other,:)];
    lonlat2=[lonlat2;pso.lonlat_other(ix_other,:)];
    clear pso slo

    pmo=load(pmothername);
    ph_patch_other2=pmo.ph_patch_other(ix_other,:);
    clear pm

    K_ps2=[K_ps2;K_ps_other2];
    C_ps2=[C_ps2;C_ps_other2];
    coh_ps2=[coh_ps2;coh_ps_other2];
    ph_patch2=[ph_patch2;ph_patch_other2];
    ph_res2=[ph_res2;ph_res_other2];
else
    n_ps_other=0;
end

if exist(hgtname,'file')
    ht=load(hgtname);
    hgt=ht.hgt(ix2);
    clear ht
    if all_da_flag~=0
        hto=load(hgtothername);
        hgt=[hgt;hto.hgt_other(ix_other)];
        clear hto
    end
end


n_ps_low_D_A=length(ix2);
n_ps=n_ps_low_D_A + n_ps_other;
ix_weed=logical(ones(n_ps,1));
disp([num2str(n_ps_low_D_A),' low D_A PS, ',num2str(n_ps_other),' high D_A PS']);




if no_weed_adjacent==0
	step_name='INITIALISE NEIGHBOUR MATRIX'
    
    ij_shift=ij2(:,2:3)+repmat([2,2]-min(ij2(:,2:3)),n_ps,1);
	neigh_ix=zeros(max(ij_shift(:,1))+1,max(ij_shift(:,2))+1);
    miss_middle=logical(ones(3));
    miss_middle(2,2)=0;
    
	for i=1:n_ps
        neigh_this=neigh_ix(ij_shift(i,1)-1:ij_shift(i,1)+1,ij_shift(i,2)-1:ij_shift(i,2)+1);
        neigh_this(neigh_this==0&miss_middle)=i;
        neigh_ix(ij_shift(i,1)-1:ij_shift(i,1)+1,ij_shift(i,2)-1:ij_shift(i,2)+1)=neigh_this;
        
        if i/100000==floor(i/100000)
            disp([num2str(i),' PS processed'])
            save log step_name i
            !sync
        end
    end
    
	step_name='FIND NEIGHBOURS'
    
    neigh_ps=cell(n_ps,1);
	for i=1:n_ps
        my_neigh_ix=neigh_ix(ij_shift(i,1),ij_shift(i,2));
        if my_neigh_ix~=0
            neigh_ps{my_neigh_ix}=[neigh_ps{my_neigh_ix},i];
        end    
        if i/100000==floor(i/100000)
            disp([num2str(i),' PS processed'])
            save log step_name i
            !sync
        end

    end	
    
    clear neigh_ix
    
	step_name='SELECT BEST'
	
	for i=1:n_ps
        if ~isempty(neigh_ps{i})
            same_ps=i;
            i2=1;
            while i2 <= length(same_ps)
                ps_i=same_ps(i2);
                same_ps=[same_ps,neigh_ps{ps_i}];
                neigh_ps{ps_i}=[];
                i2=i2+1;
            end
            same_ps=unique(same_ps);
            [dummy,high_coh]=max(coh_ps2(same_ps));
            low_coh_ix=logical(ones(size(same_ps)));
            low_coh_ix(high_coh)=0;
            ix_weed(same_ps(low_coh_ix))=0;
        end
        if i/100000==floor(i/100000)
            disp([num2str(i),' PS processed'])
            save log step_name i
            !sync
        end
	
	end
	


    disp([num2str(sum(ix_weed)),' PS kept after dropping adjacent pixels']);

end

    if strcmpi(weed_zero_elevation,'y') & exist('hgt','var')
        sea_ix=hgt==0;
        ix_weed(sea_ix)=false;
    end


    xy_weed=xy2(ix_weed,:);

    %%%%%%%%% Some non-adjacent pixels are allocated the same lon/lat by DORIS. %%%%%%%%%
    %%%%%%%%% If duplicates occur, the pixel with the highest coherence is kept.%%%%%%%%% 
    ix_weed_num=find(ix_weed); 
    [dummy,I]=unique(xy_weed(:,2:3),'rows');
    dups=setxor(I,[1:sum(ix_weed)]'); % pixels with duplicate lon/lat

    for i=1:length(dups)
        dups_ix_weed=find(xy_weed(:,2)==xy_weed(dups(i),2)&xy_weed(:,3)==xy_weed(dups(i),3));
        dups_ix=ix_weed_num(dups_ix_weed);
        [dummy,I]=max(coh_ps2(dups_ix));
        ix_weed(dups_ix([1:end]~=I))=0; % drop dups with lowest coh
    end

    if ~isempty(dups)
        xy_weed=xy2(ix_weed,:);
        fprintf('%d PS with duplicate lon/lat dropped\n\n',length(dups)')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

n_ps=sum(ix_weed);
ps_std=zeros(n_ps,1);
ps_max=zeros(n_ps,1);

if no_weed_noisy==0
    
    step_name='DROP NOISY'

    nodename=['psweed.1.node'];
    fid=fopen(nodename,'w');
    fprintf(fid,'%d 2 0 0\n',n_ps);

    for i=1:n_ps
        fprintf(fid,'%d %f %f\n',i,xy_weed(i,2),xy_weed(i,3));
    end

    fclose(fid);

    !triangle -e psweed.1.node

    fid=fopen('psweed.2.edge','r');
    header=str2num(fgetl(fid));
    N=header(1);
    edges=zeros(N,4);
    for i=1:N
        edges(i,:)=str2num(fgetl(fid));
    end
    fclose(fid);
    n_edge=size(edges,1);
    ph_weed=ph2(ix_weed,:).*exp(-j*(K_ps2(ix_weed)*bperp'));  % subtract range error 
    edge_std=zeros(n_edge,1);
    edge_max=zeros(n_edge,1);
    dph_space=(ph_weed(edges(:,3),:).*conj(ph_weed(edges(:,2),:)));

    if ~strcmpi(small_baseline_flag,'y')

      fprintf('Estimating noise for all arcs...\n')

      dph_smooth=zeros(n_edge,ps.n_ifg,'single');
      for i1=1:ps.n_ifg
        time_diff=(day(i1)-day)';
        weight_factor=exp(-(time_diff.^2)/2/time_win^2);
        weight_factor=weight_factor/sum(weight_factor);

        dph_mean=sum(dph_space.*repmat(weight_factor,n_edge,1),2);
        dph_mean_adj=angle(dph_space.*repmat(conj(dph_mean),1,ps.n_ifg)); % subtract weighted mean
        G=[ones(ps.n_ifg,1),time_diff'];
        WG=G.*[weight_factor',weight_factor'];
        m=double(WG)\(repmat(double(weight_factor),n_edge,1).*double(dph_mean_adj))'; % weighted least-sq to find best-fit line
        dph_mean_adj=angle(exp(j*(dph_mean_adj-(G*m)'))); % subtract first estimate
        m2=double(WG)\(repmat(double(weight_factor),n_edge,1).*double(dph_mean_adj))'; % weighted least-sq to find best-fit line
        dph_smooth(:,i1)=dph_mean.*exp(j*(m(1,:)'+m2(1,:)')); % add back weighted mean

       end
       dph_noise=angle(dph_space.*conj(dph_smooth));
       clear dph_space dph_smooth
       edge_std=std(dph_noise,0,2);
       edge_max=max(abs(dph_noise),[],2);
       clear dph_noise

    else
       edge_std=std(angle(dph_space),0,2);
       edge_max=max(abs(angle(dph_space)),[],2);
       clear dph_space
    end



    for i=1:n_ps
      edge_ix=[find(edges(:,2)==i);find(edges(:,3)==i)];
      ps_std(i)=min(edge_std(edge_ix)); % least noisy 
      ps_max(i)=min(edge_max(edge_ix)); % least noisy 
      if i/10000==floor(i/10000)
          fprintf('   %d PS of %d processed\n',i,n_ps)
      end
    end
    ix_weed(ix_weed)=ps_std<weed_standard_dev&ps_max<weed_max_noise;
    n_ps=sum(ix_weed);

    disp([num2str(n_ps),' PS kept after dropping noisy pixels']);

end

weedname=['weed',num2str(psver)];
save(weedname,'ix_weed','ps_std')


coh_ps=coh_ps2(ix_weed);
K_ps=K_ps2(ix_weed);
C_ps=C_ps2(ix_weed);
ph_patch=ph_patch2(ix_weed,:);
if ~isempty(ph_res2)
    ph_res=ph_res2(ix_weed,:);
else
    ph_res=ph_res2;
end

pmname=['pm',num2str(psver+1)];
save(pmname,'ph_patch','ph_res','coh_ps','K_ps','C_ps')
clear ph_patch ph_res coh_ps K_ps C_ps ph_patch2 ph_res2 coh_ps2 K_ps2 C_ps2

ph2=ph2(ix_weed,:);
ph=ph2;
phname=['ph',num2str(psver+1)];
save(phname,'ph')
clear ph

xy2=xy2(ix_weed,:);
ij2=ij2(ix_weed,:);
lonlat2=lonlat2(ix_weed,:);
ps.xy=xy2;
ps.ij=ij2;
ps.lonlat=lonlat2;
ps.n_ps=size(ph2,1);
psname=['ps',num2str(psver+1)];
save(psname,'-struct','ps');
clear ps xy2 ij2 lonlat2

if exist(hgtname,'file')
    hgt=hgt(ix_weed);
    save(['hgt',num2str(psver+1),'.mat'],'hgt');
    clear hgt
end

if exist(laname,'file')
    la=load(laname);
    la=[la.la(ix2)];
    if all_da_flag~=0
        lao=load(laothername);
        la=[la;lao.la_other(ix_other)];
        clear lao
    end
    la=la(ix_weed);
    save(['la',num2str(psver+1),'.mat'],'la');
    clear la
end

if exist(bpname,'file')
    bp=load(bpname);
    bperp_mat=[bp.bperp_mat(ix2,:)];
    clear bp
    if all_da_flag~=0
        bpo=load(bpothername);
        bperp_mat=[bperp_mat;bpo.bperp_other(ix_other,:)];
        clear bpo
    end
    bperp_mat=bperp_mat(ix_weed,:);
    save(['bp',num2str(psver+1),'.mat'],'bperp_mat');
end

if exist(['scla_smooth',num2str(psver+1),'.mat'],'file')
    delete(['scla_smooth',num2str(psver+1),'.mat'])
end
if exist(['scla',num2str(psver+1),'.mat'],'file')
    delete(['scla',num2str(psver+1),'.mat'])
end
if exist(['scla_smooth_sb',num2str(psver+1),'.mat'],'file')
    delete(['scla_smooth_sb',num2str(psver+1),'.mat'])
end
if exist(['scla_sb',num2str(psver+1),'.mat'],'file')
    delete(['scla_sb',num2str(psver+1),'.mat'])
end
if exist(['aps',num2str(psver+1),'.mat'],'file')
    delete(['aps',num2str(psver+1),'.mat'])
end
if exist(['scn',num2str(psver+1),'.mat'],'file')
    delete(['scn',num2str(psver+1),'.mat'])
end

setpsver(psver+1)
