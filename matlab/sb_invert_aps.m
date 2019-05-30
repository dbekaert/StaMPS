function []=sb_invert_aps(aps_flag)
%SB_invert_aps Invert unwrapped phase of short baseline ifgs
%
%   Bekaert David, November 2013   - University of Leeds
%
% modifications:
% 11/2013   DB  Include WRF model and change dry to hydrostatic
% 03/2014   AH  Add backwards compatibility 
% 12/2014   DB  Invert such that when one date exist the delay is excluded
% 09/2015   DB  Add some extra options such as MODIS. 
% 08/2017   DB  Using aps_save to store the variables
% 10/2017   DB  Bug fix for error catching of missing master
% 11/2017   DB  Adding merra, merra2, and gacos model
% 02/2018   DB  Adding in the NARR model

logit;

load psver
psname=['./ps',num2str(psver)];
apssbname=['./tca_sb',num2str(psver) '.mat'];
apsname=['./tca',num2str(psver) '.mat'];


if nargin<1
   aps_flag = []; 
end
ps=load(psname);

drop_ifg_index=getparm('drop_ifg_index');
unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);

ref_ps=ps_setref;

G=zeros(ps.n_ifg,ps.n_image);
for i=1:ps.n_ifg
    G(i,ps.ifgday_ix(i,1))=-1;
    G(i,ps.ifgday_ix(i,2))=1;
end
if sum(abs(G(:,ps.master_ix)))==0; 
   error('Apparently none of the unwrapped interferograms include the original master image')
else
    G(:,ps.master_ix)=0; % take out master as ref by setting to zero
end


G2=G(unwrap_ifg_index,:);
nzc_ix=sum(abs(G2))~=0; % index for non-zero columns
G2=G2(:,nzc_ix);
if rank(G2)<size(G2,2) 
    error('There are isolated subsets (cannot be inverted w.r.t. master)')
end



if ~isempty(aps_flag)
    aps=load(apssbname);
    [aps_corr_sb,fig_name_tca,aps_flag] = ps_plot_tca(aps,aps_flag);
    ix = find(nanmean(aps_corr_sb,1)==0);
    unwrap_ifg_index_new = setdiff(unwrap_ifg_index,ix);

    G3=G(unwrap_ifg_index_new,:);
    nzc_ix=sum(abs(G3))~=0; % index for non-zero columns
    G3=G3(:,nzc_ix);
    if rank(G3)<size(G3,2) 
        error('There are isolated subsets that do not have an APS estimate (cannot be inverted w.r.t. master)')
    end

    aps_corr_sb= aps_corr_sb(:,unwrap_ifg_index_new);
    G2 =G3;
    
    % check if the master to be inverted actually has a delay estimated with it.
    ifgday_new = ps.ifgday(unwrap_ifg_index_new,:);
    ifgs_days = unique(reshape(ifgday_new((sum(aps_corr_sb)~=0),:),[],1));

    if sum(ps.master_day==ifgs_days)~=1
        error('Master does not have an APS estimated, inversion not possible')
    end
    
    aps_corr=zeros(ps.n_ps,ps.n_image,'single');
    aps_corr(:,nzc_ix)=lscov(G2,double(aps_corr_sb'))';
    
    
    
    
    if aps_flag==1 % linear correction
        ph_tropo_linear = aps_corr;       
        aps_save(apsname,ph_tropo_linear)
    elseif aps_flag==2 % powerlaw correlation
        ph_tropo_powerlaw = aps_corr;       
        aps_save(apsname,ph_tropo_powerlaw)
    elseif aps_flag==3 % meris correction
        ph_tropo_meris = aps_corr;       
        aps_save(apsname,ph_tropo_meris)        
    elseif aps_flag==4 % ERA-I correction
        ph_tropo_era = aps_corr;
        aps_save(apsname,ph_tropo_era)       
    elseif aps_flag==5 % ERA-I correction
        ph_tropo_era_hydro = aps_corr;
        aps_save(apsname,ph_tropo_era_hydro)       
    elseif aps_flag==6 % ERA-I correction
        ph_tropo_era_wet = aps_corr;
        aps_save(apsname,ph_tropo_era_wet)       
    elseif aps_flag==7 % WRF correction
        ph_tropo_wrf = aps_corr;
        aps_save(apsname,ph_tropo_wrf)       
    elseif aps_flag==8 % WRF correction
        ph_tropo_wrf_hydro = aps_corr;
        aps_save(apsname,ph_tropo_wrf_hydro)       
    elseif aps_flag==9 % WRF correction
        ph_tropo_wrf_wet = aps_corr;
        aps_save(apsname,ph_tropo_wrf_wet)       
    elseif aps_flag==12 % modis correction
        ph_tropo_modis= aps_corr;
        aps_save(apsname,ph_tropo_modis)

    elseif aps_flag==13 % modis correction (not interpolated)
        ph_tropo_modis_no_interp= aps_corr;
        aps_save(apsname,ph_tropo_modis_no_interp)       

     elseif aps_flag==18 % current implementation of aps correction (manually estimated)
        strat_corr = aps_corr; % Old implementation
        aps_save(apsname,strat_corr)
      
    elseif aps_flag==19 % modis correction
        ph_tropo_modis_recal= aps_corr;
        aps_save(apsname,ph_tropo_modis_recal)
    elseif aps_flag==20 % modis correction (not interpolated)
        ph_tropo_modis_no_interp_recal= aps_corr;
        aps_save(apsname,ph_tropo_modis_no_interp_recal)   
        
        
        
     elseif aps_flag==29 % MERRA correction
        ph_tropo_merra = aps_corr;
        aps_save(apsname,ph_tropo_merra)       
    elseif aps_flag==31 % MERRA correction
        ph_tropo_merra_hydro = aps_corr;
        aps_save(apsname,ph_tropo_merra_hydro)       
    elseif aps_flag==33 % MERRA correction
        ph_tropo_merra_wet = aps_corr;
        aps_save(apsname,ph_tropo_merra_wet)     
        
    elseif aps_flag==30 % MERRA-2 correction
        ph_tropo_merra2 = aps_corr;
        aps_save(apsname,ph_tropo_merra2)       
    elseif aps_flag==32 % MERRA-2 correction
        ph_tropo_merra2_hydro = aps_corr;
        aps_save(apsname,ph_tropo_merra2_hydro)       
    elseif aps_flag==34 % MERRA-2 correction
        ph_tropo_merra2_wet = aps_corr;
        aps_save(apsname,ph_tropo_merra2_wet)    
        
   elseif aps_flag==35 % GACOS correction
        ph_tropo_gacos = aps_corr;
        aps_save(apsname,ph_tropo_gacos)        
   elseif aps_flag==36 % NARR correction
        ph_tropo_narr = aps_corr;
        aps_save(apsname,ph_tropo_narr)   
   elseif aps_flag==37 % NARR correction hydro
        ph_tropo_narr_hydro = aps_corr;
        aps_save(apsname,ph_tropo_narr_hydro)   
   elseif aps_flag==38 % NARR correction wet
        ph_tropo_narr_wet = aps_corr;
        aps_save(apsname,ph_tropo_narr_wet) 
   elseif aps_flag==39 % ERA-5 correction
        ph_tropo_era5 = aps_corr;
        aps_save(apsname,ph_tropo_era5)       
   elseif aps_flag==40 % ERA-5 correction
        ph_tropo_era5_hydro = aps_corr;
        aps_save(apsname,ph_tropo_era5_hydro)       
   elseif aps_flag==41 % ERA-5 correction
        ph_tropo_era5_wet = aps_corr;
        aps_save(apsname,ph_tropo_era5_wet) 
   else
        error('Currently other options not supported.\nIf you are combining a combination of two different techniques, try the following:\n plot each technique for a SM correction\n')
    end
end



logit(1);
