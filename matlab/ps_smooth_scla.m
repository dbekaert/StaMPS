function []=ps_smooth_scla(use_small_baselines)
%PS_SMOOTH_SCLA 
%
%   Andy Hooper, March 2007
%
%   ======================================================================
%   03/2009 AH: save in scla_smooth mat files
%   06/2009 AH: orbital ramps added 
%   ======================================================================
%
fprintf('Smoothing spatially-correlated look angle error...\n')

if nargin<1
    use_small_baselines=0;
end

scn_wavelength=100;
small_baseline_flag=getparm('small_baseline_flag');

load psver
psname=['ps',num2str(psver)];
bpname=['bp',num2str(psver)];
if use_small_baselines==0
    sclaname=['scla',num2str(psver)];
    sclasmoothname=['scla_smooth',num2str(psver)];
else
    sclaname=['scla_sb',num2str(psver)];
    sclasmoothname=['scla_smooth_sb',num2str(psver)];
end

ps=load(psname);
scla=load(sclaname);

n_ps=ps.n_ps;
n_ifg=ps.n_ifg;

disp(sprintf('   Number of points per ifg: %d',n_ps))

sigma_sq_times_2=2*scn_wavelength.^2;
patch_dist=scn_wavelength*4;
patch_dist_sq=patch_dist*patch_dist;
ix_range=ceil(n_ps/(max(ps.xy(:,3))-min(ps.xy(:,3)))*patch_dist*0.2);
ix1=1;
ix2=ix_range;
ps.xy(:,1)=[1:n_ps]';
K_ps_uw_smooth=zeros(size(scla.K_ps_uw));
C_ps_uw_smooth=K_ps_uw_smooth;

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
    K_ps_uw_smooth(i)=sum(scla.K_ps_uw(xy_near(:,1)).*weight_factor);
    C_ps_uw_smooth(i)=sum(scla.C_ps_uw(xy_near(:,1)).*weight_factor);

    if i/10000==floor(i/10000)
        fprintf('   %d PS processed\n',i)
    end % end-if
end % i++

K_ps_uw=K_ps_uw_smooth;
C_ps_uw=C_ps_uw_smooth;
ph_ramp=scla.ph_ramp;

clear scla

bp=load(bpname);
if use_small_baselines==0
    if strcmpi(small_baseline_flag,'y')
        bperp_mat=zeros(ps.n_ps,ps.n_image-1);
        G=zeros(ps.n_ifg,ps.n_image);
        for i=1:ps.n_ifg
             G(i,ps.ifgday_ix(i,1))=-1;
             G(i,ps.ifgday_ix(i,2))=1;
        end
        G=G(:,[1:ps.master_ix-1,ps.master_ix+1:end]);
        bperp_mat=[G\double(bp.bperp_mat')]';
        bperp_mat=[bperp_mat(:,1:ps.master_ix-1),zeros(ps.n_ps,1,'single'),bperp_mat(:,ps.master_ix:end)];
        %unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
    else

        bperp_mat=[bp.bperp_mat(:,1:ps.master_ix-1),zeros(ps.n_ps,1,'single'),bp.bperp_mat(:,ps.master_ix:end)];
        %unwrap_ifg_index=setdiff(unwrap_ifg_index,ps.master_ix);
    end
else
    bperp_mat=bp.bperp_mat;
end

ph_scla=repmat(K_ps_uw,1,size(bperp_mat,2)).*bperp_mat;

save(sclasmoothname,'K_ps_uw','C_ps_uw','ph_scla','ph_ramp')    
