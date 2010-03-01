function []=ps_scn_filt()
%PS_SCN_FILT estimate spatially correlated noise in unwrapped phase
%
%   Andy Hooper, June 2006
%
%   =======================================================================
%   11/2006 AH: Error corrected that was leaving master in temporal smoothing
%   04/2007 AH: Added 64-bit machine compatibility
%   05/2007 AH: Spatially correlated look angle error added  
%   =======================================================================

fprintf('Estimating other spatially-correlated noise...\n')

pix_size=getparm('unwrap_grid_size',1);
time_win=getparm('scn_time_win',1);
deramp_ifg=getparm('scn_deramp_ifg',1);
scn_wavelength=getparm('scn_wavelength',1);
unwrap_ifg_index=getparm('unwrap_ifg_index',1);
small_baseline_flag=getparm('small_baseline_flag',1);

load psver
psname=['ps',num2str(psver)];
phuwname=['phuw',num2str(psver)];
sclaname=['scla',num2str(psver)];
apsname=['aps',num2str(psver)];
scnname=['scn',num2str(psver)]; % spatially-correlated noise

ps=load(psname);
uw=load(phuwname);
%aps=load(apsname);

if strcmpi(small_baseline_flag,'y')
    unwrap_ifg_index=[1:ps.n_image];
elseif strcmp(unwrap_ifg_index,'all')
    unwrap_ifg_index=[1:ps.n_ifg];
end

day=ps.day(unwrap_ifg_index);
master_ix=sum(ps.master_day>day)+1;
n_ifg=length(unwrap_ifg_index);
n_ps=ps.n_ps;

ph_all=single(uw.ph_uw(:,unwrap_ifg_index));
if exist([sclaname,'.mat'],'file')
    scla=load(sclaname);
    ph_all=ph_all-single(scla.ph_scla(:,unwrap_ifg_index));
    ph_all=ph_all-repmat(single(scla.C_ps_uw),1,length(unwrap_ifg_index));
end
ph_all(isnan(ph_all))=0;  

disp(sprintf('   Number of points per ifg: %d',n_ps))

nodename=['scnfilt.1.node'];
fid=fopen(nodename,'w');
fprintf(fid,'%d 2 0 0\n',n_ps);

for i=1:n_ps
    fprintf(fid,'%d %f %f\n',i,ps.xy(i,2),ps.xy(i,3));
end

fclose(fid);

!triangle -e scnfilt.1.node

fid=fopen('scnfilt.2.edge','r');
header=str2num(fgetl(fid));
N=header(1);
edges_nz=zeros(N,4);
for i=1:N
    edges_nz(i,:)=str2num(fgetl(fid));
end
fclose(fid);

%%% deramp end ifgs (unlike aps, orbit errors not so random and end
%%% orbit errors can pass through the low-pass filter
deramp_ifg=intersect(deramp_ifg,unwrap_ifg_index);
deramp_ix=zeros(size(deramp_ifg));
ph_ramp=zeros(n_ps,length(deramp_ifg));

if ~isempty(deramp_ifg)
    fprintf('   deramping selected ifgs...\n')
    G=double([ones(n_ps,1),ps.xy(:,2),ps.xy(:,3)]);
    %G=double([ones(n_ps,1),ps.xy(:,2)]); % range only

    for i=1:length(deramp_ifg)
        i3=find(unwrap_ifg_index==deramp_ifg(i))
        deramp_ix(i)=i3;
        d=(ph_all(:,i3));
        m=G\double(d(:));
        ph_this_ramp=G*m;
        ph_all(:,i3)=ph_all(:,i3)-ph_this_ramp; % subtract ramp
        ph_ramp(:,i)=ph_this_ramp;
    end
    save(scnname,'ph_ramp')    
end


%%% smooth in time using gaussian moving window
isnanix=isnan(uw.ph_uw);
uw.ph_uw(isnanix)=0;
dph=ph_all(edges_nz(:,3),:)-ph_all(edges_nz(:,2),:);
dph_lpt=zeros(size(dph));
n_edges=size(dph,1);

fprintf('   low-pass filtering pixel-pairs in time...\n')

for i1=1:n_ifg
    time_diff_sq=(day(i1)-day)'.^2;
    weight_factor=exp(-time_diff_sq/2/time_win^2);
    weight_factor(master_ix)=0; % leave out master
    weight_factor=weight_factor/sum(weight_factor);
    dph_lpt(:,i1)=sum(dph.*repmat(weight_factor,n_edges,1),2);
end


dph_hpt=dph-dph_lpt;  % leaves master APS - slave APS - slave noise (+ residue master noise)

ph_hpt=zeros(n_ps-1,n_ifg);
ref_ix=1;

A=sparse([[1:n_edges]';[1:n_edges]'],[edges_nz(:,2);edges_nz(:,3)],[-ones(n_edges,1);ones(n_edges,1)]);
A=double(A(:,[1:ref_ix-1,ref_ix+1:n_ps]));

fprintf('   solving for high-frequency (in time) pixel phase...\n')

for i=1:n_ifg
    ph_hpt(:,i)=A\double(dph_hpt(:,i));
end

ph_hpt=[ph_hpt(1:ref_ix-1,:);zeros(1,n_ifg);ph_hpt(ref_ix:end,:)]; % add back ref point


ph_hpt(:,deramp_ix)=ph_hpt(:,deramp_ix)+ph_ramp;

ph_hpt=single(ph_hpt);

sigma_sq_times_2=2*scn_wavelength.^2;
ph_scn=nan(n_ps,n_ifg);
patch_dist=scn_wavelength*4;
patch_dist_sq=patch_dist*patch_dist;
ix_range=ceil(n_ps/(max(ps.xy(:,3))-min(ps.xy(:,3)))*patch_dist*0.2);
ix1=1;
ix2=ix_range;
ps.xy(:,1)=[1:n_ps]';

fprintf('   low-pass filtering in space...\n')

for i=1:n_ps
    x_min=ps.xy(i,2)-patch_dist;
    x_max=ps.xy(i,2)+patch_dist;
    y_min=ps.xy(i,3)-patch_dist;
    y_max=ps.xy(i,3)+patch_dist;

    ix1=ix1+ix_range;
    ix1(ix1>n_ps)=n_ps;
    while ix1>1 & ps.xy(ix1-1,3)>=y_min
        ix1=ix1-ix_range;
    end

    ix2=ix2-ix_range;
    ix2(ix2<1)=1;
        while ix2<n_ps & ps.xy(ix2+1,3)<=y_max
        ix2=ix2+ix_range;
    end

    ix1(ix1<1)=1;
    ix2(ix2>n_ps)=n_ps;

    xy_near=ps.xy(ix1:ix2,:);
    xy_near=xy_near(xy_near(:,2)>=x_min & xy_near(:,2)<=x_max & xy_near(:,3)>=y_min & xy_near(:,3)<=y_max,:);
    dist_sq=(xy_near(:,2)-ps.xy(i,2)).^2+(xy_near(:,3)-ps.xy(i,3)).^2;
    in_range_ix=dist_sq<patch_dist_sq; % exclude those out of range
    xy_near=xy_near(in_range_ix);
    dist_sq=dist_sq(in_range_ix); 
    weight_factor=exp(-dist_sq/sigma_sq_times_2);
    weight_factor=weight_factor/sum(weight_factor); % normalize
    ph_scn(i,:)=sum(ph_hpt(xy_near(:,1),:).*repmat(weight_factor,1,n_ifg),1);

    if i/1000==floor(i/1000)
        disp([num2str(i),' PS processed'])
    end % end-if
end % i++

ph_scn=ph_scn-repmat(ph_scn(1,:),n_ps,1); % re-ref to 1st PS
ph_scn_slave=zeros(size(uw.ph_uw));
ph_scn_slave(:,unwrap_ifg_index)=ph_scn;
ph_scn_slave(:,master_ix)=0;

save(scnname,'ph_scn_slave','ph_hpt','ph_ramp')    
