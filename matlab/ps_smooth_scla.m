function []=ps_smooth_scla(use_small_baselines)
%PS_SMOOTH_SCLA 
%
%   Andy Hooper, March 2007
%
%   ======================================================================
%   03/2009 AH: save in scla_smooth mat files
%   06/2009 AH: orbital ramps added 
%   01/2012 AH: Filtering strategy changed to just remove outliers
%   09/2015 AH: use matlab triangulation if triangle program not installed
%   06/2017 DB: Include stamps save for large variables
%   ======================================================================
%
logit;
logit('Smoothing spatially-correlated look angle error...',2)

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
K_ps_uw=scla.K_ps_uw;
C_ps_uw=scla.C_ps_uw;
ph_ramp=scla.ph_ramp;
clear scla

n_ps=ps.n_ps;
n_ifg=ps.n_ifg;

logit(sprintf('Number of points per ifg: %d',n_ps))

arch=computer('arch');
if strcmpi(arch(1:3),'win')
    use_triangle='n';
else
    tripath=system('which triangle >& /dev/null');
    if tripath==0
        use_triangle='y';
    else
        use_triangle='n';
    end  
end
    
if use_triangle=='y'
    nodename=['scla.1.node'];
    fid=fopen(nodename,'w');
    fprintf(fid,'%d 2 0 0\n',n_ps);
    ps.xy(:,1)=1:n_ps;
    fprintf(fid,'%d %f %f\n',ps.xy');
    fclose(fid);

    system('triangle -e scla.1.node > triangle_scla.log');

    fid=fopen('scla.2.edge','r');
    header=str2num(fgetl(fid));
    N=header(1);
    edgs=zeros(N,4);
    edgs=fscanf(fid,'%d %d %d %d\n',[4,N])';
    fclose(fid);
    edgs=edgs(:,2:3);
    n_edge=size(edgs,1);
    if n_edge~=N
        error('missing lines in scla.2.edge')
    end
else
    xy=double(ps.xy);
    tri=delaunay(xy(:,2),xy(:,3));
    tr=triangulation(tri,xy(:,2),xy(:,3));
    edgs=edges(tr);
    n_edge=size(edgs,1);
end

  
logit(sprintf('Number of arcs per ifg=%d',n_edge))

Kneigh_min=inf(n_ps,1,'single');
Kneigh_max=-inf(n_ps,1,'single');
Cneigh_min=inf(n_ps,1,'single');
Cneigh_max=-inf(n_ps,1,'single');

for i=1:n_edge % find min and max neighbour for each ps
    ix=edgs(i,1:2);
    Kneigh_min(ix)=min([Kneigh_min(ix),K_ps_uw(fliplr(ix))],[],2);
    Kneigh_max(ix)=max([Kneigh_max(ix),K_ps_uw(fliplr(ix))],[],2);
    Cneigh_min(ix)=min([Cneigh_min(ix),C_ps_uw(fliplr(ix))],[],2);
    Cneigh_max(ix)=max([Cneigh_max(ix),C_ps_uw(fliplr(ix))],[],2);
    if i/100000==floor(i/100000)
        logit(sprintf('%d arcs processed',i),2)
    end % end-if
end
    
ix1=K_ps_uw>Kneigh_max;
ix2=K_ps_uw<Kneigh_min;
K_ps_uw(ix1)=Kneigh_max(ix1); % reduce positive outliers
K_ps_uw(ix2)=Kneigh_min(ix2); % increase negative outliers
   
ix1=C_ps_uw>Cneigh_max;
ix2=C_ps_uw<Cneigh_min;
C_ps_uw(ix1)=Cneigh_max(ix1); % reduce positive outliers
C_ps_uw(ix2)=Cneigh_min(ix2); % increase negative outliers

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
    else
        bperp_mat=[bp.bperp_mat(:,1:ps.master_ix-1),zeros(ps.n_ps,1,'single'),bp.bperp_mat(:,ps.master_ix:end)];
    end
else
    bperp_mat=bp.bperp_mat;
end

ph_scla=repmat(K_ps_uw,1,size(bperp_mat,2)).*bperp_mat;

stamps_save(sclasmoothname,K_ps_uw,C_ps_uw,ph_scla,ph_ramp)    
logit(1);
