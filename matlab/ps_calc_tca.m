function []=ps_calc_tca
%PS_CALC_TCA calculate topo-correlated atmosphere
%   execute in according INSAR (for PS) or SMALL_BASELINES (for SB) ifg
%   correction after running stamps(1,7)
%
%   The script will load already existing correction 
%   Keep in mind that it is better to apply no
%   correction if the correlation is low !
%
%   Hannes Bathke, November 2012
%
%   ================================================================
%   11/2012 HB added more precise pixel selection, bugfixes
%   11/2012 HB drop ifg index included
%   01/2013 AH converted to function and name changed
%   08/2014 DB Append the correction in case the tca variable esits 
%   08/2014 DB Change such can be ran without having ran step 7
%   ================================================================

clear all
close all

plot_flag = 1;    % 0 --> no plot/ 1 -> plot uph, correction, corrected ifg after each correction
%% load standard ps indata check if sbas or ps folder
% load interferograms

load psver
psname=['./ps',num2str(psver)];
sclaname=['./scla',num2str(psver) '.mat'];
sclasbname=['./scla_sb',num2str(psver) '.mat'];
phuwname=['./phuw',num2str(psver) '.mat'];
phuwsbname=['./phuw_sb',num2str(psver) '.mat'];


ps=load(psname);

num_ifg = ps.n_ifg;
num_pix = ps.n_ps;

small_baseline_flag = getparm('small_baseline_flag');
if ~strcmpi(small_baseline_flag,'y')
    disp('working on the PS interferograms');
    if exist([phuwname],'file')
        stratname=['./tca',num2str(psver)]; 
        u_ph = load(phuwname');
        
    end
    if exist(sclasbname,'file')==2
        scla=load(sclaname);
        % subtract dem error
        ph_all=u_ph.ph_uw - scla.ph_scla;
        
        fprintf('DEM error removed \n')
    else
        ph_all=u_ph.ph_uw  ;      
    end
else
    disp('working on the SB interferograms');
    if exist([phuwsbname],'file')
        stratname=['./tca_sb',num2str(psver)]; 
        u_ph = load(phuwsbname);  
     end
    if exist(sclasbname,'file')==2
        scla=load(sclasbname);
        % subtract dem error
        ph_all=u_ph.ph_uw - scla.ph_scla;
        fprintf('DEM error removed \n')

    else
        ph_all=u_ph.ph_uw;
    end
end

% if existent load strat correction
strat_corr = zeros(num_pix,num_ifg);


% load dem data
if exist('hgt2.mat','file')
    dem = load('hgt2.mat');
else
    disp('dem height not existing yet, please run StaMPS(1,7) once')
    exit
end

% select interferograms for correction
drop_ifg_index = getparm('drop_ifg_index');
disp(['dropped ifgs: ',num2str(drop_ifg_index)]);
disp(['Loop through ',num2str(num_ifg),' ifgs:']);
ifg = input('Start with which ifg? ');
end_ifg = input('End with which ifg? ');

for k = ifg : end_ifg
    disp(['interferogram number ',num2str(k)])
    % check if ifg is on drop index
    if k == drop_ifg_index;
        disp('ifg dropped, nothing done, next');
    else
        figure;
        plot(dem.hgt,ph_all(:,k),'.k')
        hold on
        xlabel('DEM Elevation [m]');
        ylabel('unwrapped phase')
        do_corr = input('correlation existent? (yes - 1; no(next ifg) - 0 ) ');
        if do_corr == 1
            disp('Click 2 edges of a rectangle within the linear function\n should be fitted to the pixels ');
            bounds = ginput(2);
            low_hgt = min(bounds(:,1));
            high_hgt = max(bounds(:,1));
            low_uph = min(bounds(:,2));
            high_uph = max(bounds(:,2));
            %find pixels within bounds
            did = find(low_hgt < dem.hgt & dem.hgt < high_hgt);
            phid = find(low_uph < ph_all(:,k) & ph_all(:,k) < high_uph);
            pix_int = intersect(did,phid);
            plot(dem.hgt(pix_int),ph_all(pix_int,k),'.b')
            
            % build design matrix
            dem_range = low_hgt:1:high_hgt;
            A = [dem.hgt(pix_int) ones(length(pix_int),1)];
            Afull = [dem_range' ones(length(dem_range),1)];
            
            % do the inversion
            coef = (A'*A)\A'*ph_all(pix_int,k);
            strat_fun = Afull * coef;
            
            % display fitted function
            plot(dem_range,strat_fun,'-r')
            hold off
            strat_corr(:,k) = coef(1)*dem.hgt;
            
            if (plot_flag)
                % plot results
                figure; scatter(ps.lonlat(:,1),ps.lonlat(:,2),5,strat_corr(:,k),'filled');
                title('estimated stratified atm comp'); colorbar
                caxis([min(ph_all(:,k)) max(ph_all(:,k))])
                
                figure; scatter(ps.lonlat(:,1),ps.lonlat(:,2),5,ph_all(:,k),'filled');
                title(['unwrapped phase of ifg ',num2str(k)]); colorbar
                caxis([min(ph_all(:,k)) max(ph_all(:,k))])
                
                figure; scatter(ps.lonlat(:,1),ps.lonlat(:,2),5,...
                    (ph_all(:,k)-strat_corr(:,k)),'filled');
                title(['corrected ifg ',num2str(k)]); colorbar
                caxis([min(ph_all(:,k)) max(ph_all(:,k))])
                
                next = input('next ifg? (1 --> yes) (0 --> EXIT) ');
                if next == 1
                    close all
                else
                    break
                end
            end
            close all
        else
            disp(['no correlation in interferogram ',num2str(k),' no correction applied--> 0']);
            close all
        end
    end
end

% save matrix
if exist([stratname,'.mat'],'file')~=2
    save([stratname,'.mat'],'strat_corr')
else
    save([stratname,'.mat'],'-append','strat_corr')
end
