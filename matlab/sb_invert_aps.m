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
    ifgs_days = unique(reshape(ps.ifgday((sum(aps_corr_sb)~=0),:),[],1));
    if sum(ps.master_day==ifgs_days)~=1
%         fprintf('************ \nMaster day is not esimated in the APS \n The interferometric delays will actually be SAR delays! \n *************** \n')
        error('Master does not have an APS estimated, inversion not possible')
    end
    
    aps_corr=zeros(ps.n_ps,ps.n_image,'single');
    aps_corr(:,nzc_ix)=lscov(G2,double(aps_corr_sb'))';
    
    
    
    
    if aps_flag==1 % linear correction
        ph_tropo_linear = aps_corr;       
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_linear')
        else
            save(apsname,'ph_tropo_linear')
        end  
    elseif aps_flag==2 % powerlaw correlation
        ph_tropo_powerlaw = aps_corr;       
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_powerlaw')
        else
            save(apsname,'ph_tropo_powerlaw')
        end  
    elseif aps_flag==3 % meris correction
         ph_tropo_meris = aps_corr;       
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_meris')
        else
            save(apsname,'ph_tropo_meris')        
        end 
    elseif aps_flag==4 % ERA-I correction
        ph_tropo_era = aps_corr;
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_era')
        else
            save(apsname,'ph_tropo_era')       
        end
    elseif aps_flag==5 % ERA-I correction
        ph_tropo_era_hydro = aps_corr;
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_era_hydro')
        else
            save(apsname,'ph_tropo_era_hydro')       
        end
    elseif aps_flag==6 % ERA-I correction
        ph_tropo_era_wet = aps_corr;
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_era_wet')
        else
            save(apsname,'ph_tropo_era_wet')       
        end
    elseif aps_flag==7 % WRF correction
        ph_tropo_wrf = aps_corr;
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_wrf')
        else
            save(apsname,'ph_tropo_wrf')       
        end
    elseif aps_flag==8 % WRF correction
        ph_tropo_wrf_hydro = aps_corr;
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_wrf_hydro')
        else
            save(apsname,'ph_tropo_wrf_hydro')       
        end
    elseif aps_flag==9 % WRF correction
        ph_tropo_wrf_wet = aps_corr;
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_wrf_wet')
        else
            save(apsname,'ph_tropo_wrf_wet')       
        end
        
        
    elseif aps_flag==12 % modis correction
        ph_tropo_modis= aps_corr;
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_modis')
        else
            save(apsname,'ph_tropo_modis')       
        end
    elseif aps_flag==13 % modis correction (not interpolated)
        ph_tropo_modis_no_interp= aps_corr;
         if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_modis_no_interp')
        else
            save(apsname,'ph_tropo_modis_no_interp')       
         end

      
        
     elseif aps_flag==18 % current implementation of aps correction (manually estimated)
        strat_corr = aps_corr; % Old implementation
        if exist(apsname,'file')==2
            save(apsname,'-append','strat_corr')
        else
            save(apsname,'strat_corr')       
        end
        
    elseif aps_flag==19 % modis correction
        ph_tropo_modis_recal= aps_corr;
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_modis_recal')
        else
            save(apsname,'ph_tropo_modis_recal')       
        end

    elseif aps_flag==20 % modis correction (not interpolated)
        ph_tropo_modis_no_interp_recal= aps_corr;
        if exist(apsname,'file')==2
            save(apsname,'-append','ph_tropo_modis_no_interp_recal')
        else
            save(apsname,'ph_tropo_modis_no_interp_recal')       
        end
    else
        error('Currently other options not supported.\nIf you are combining a combination of two different techniques, try the following:\n plot each technique for a SM correction\n')
    end
end



logit(1);
