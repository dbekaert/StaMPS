function []=ps_output()
%PS_OUTPUT write various output files 
%
%   Andy Hooper, June 2006
%
%   =======================================================================
%   09/2009 AH: Correct processing for small baselines output
%   =======================================================================

fprintf('Writing output files...\n')

small_baseline_flag=getparm('small_baseline_flag',1);

load psver
psname=['ps',num2str(psver)];
rcname=['rc',num2str(psver)];
phuwname=['phuw',num2str(psver)];
sclaname=['scla',num2str(psver)];
hgtname=['hgt',num2str(psver)];
scnname=['scn',num2str(psver)];
meanvname=['mean_v'];

ps=load(psname);
phuw=load(phuwname);
rc=load(rcname);

if strcmpi(small_baseline_flag,'y')
    n_image=ps.n_image;
else
    n_image=ps.n_ifg;
end

ijname=['ps_ij.txt'];
ij=ps.ij(:,2:3);
save(ijname,'ij','-ASCII');

llname=['ps_ll.txt'];
lonlat=ps.lonlat;
save(llname,'lonlat','-ASCII');

datename=['date.txt'];
date_out=str2num(datestr(ps.day,'yyyymmdd'));
save(datename,'date_out','-ASCII','-DOUBLE');

master_ix=sum(ps.master_day>ps.day)+1;

ref_ps=ps_setref;
ph_uw=phuw.ph_uw-repmat(mean(phuw.ph_uw(ref_ps,:)),ps.n_ps,1);
ph_w=angle(rc.ph_rc.*repmat(conj(sum(rc.ph_rc(ref_ps,:))),ps.n_ps,1));
ph_w(:,master_ix)=0;


fid=fopen('ph_w.flt','w');
fwrite(fid,ph_w','float');
fclose(fid);

fid=fopen('ph_uw.flt','w');
fwrite(fid,ph_uw','float');
fclose(fid);

scla=load(sclaname);
if exist(hgtname,'file')
    hgt=load(hgtname);
else
    hgt.hgt=zeros(ps.n_ps,1);
end

ph_uw=phuw.ph_uw - scla.ph_scla - repmat(scla.C_ps_uw,1,n_image);

%%% this is only approximate
K_ps_uw=scla.K_ps_uw-mean(scla.K_ps_uw);
dem_error=double(K2q(K_ps_uw,ps.ij(:,3)));
dem_error=dem_error-mean(dem_error(hgt.hgt==0));
dem_sort=sort(dem_error);
min_dem=dem_sort(ceil(length(dem_sort)*0.001));
max_dem=dem_sort(floor(length(dem_sort)*0.999));
dem_error_tt=dem_error;
dem_error_tt(dem_error<min_dem)=min_dem; % for plotting purposes
dem_error_tt(dem_error>max_dem)=max_dem; % for plotting purposes
dem_error_tt=[ps.lonlat,dem_error_tt];

save('dem_error.xy','dem_error_tt','-ASCII');

%%%

clear scla phuw
ph_uw=ph_uw-repmat(mean(ph_uw(ref_ps,:)),ps.n_ps,1);

meanv=load(meanvname);
lambda=getparm('lambda');
mean_v=-meanv.m(2,:)'*365.25/4/pi*lambda*1000; % m(1,:) is master APS + mean deviation from model
v_sort=sort(mean_v);
min_v=v_sort(ceil(length(v_sort)*0.001))
max_v=v_sort(floor(length(v_sort)*0.999))
mean_v(mean_v<min_v)=min_v;
mean_v(mean_v>max_v)=max_v;


mean_v_name=['ps_mean_v.xy'];
mean_v=[ps.lonlat,double(mean_v)];
save(mean_v_name,'mean_v','-ASCII');

%%Note mean_v is relative to a reference point
%%and dem_error is relative to mean of zero height points (if there are any)
fid=fopen('ps_data.xy','w');
fprintf(fid,'%f %f %4.4f %4.4f %4.4f\n',[mean_v,double(hgt.hgt),dem_error]');
fclose(fid)

for i=1:n_image
    ph=ph_uw(:,i);
    ph_sort=sort(ph);
    min_ph=ph_sort(ceil(length(ph_sort)*0.001));
    max_ph=ph_sort(floor(length(ph_sort)*0.999));
    ph(ph<min_ph)=min_ph;
    ph(ph>max_ph)=max_ph;
    ph=-ph*lambda*1000/4/pi;
    ph=[ps.lonlat,double(ph)];
    save(['ps_u-dm.',num2str(i),'.xy'],'ph','-ASCII');
end

