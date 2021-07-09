function stamps(start_step,end_step,patches_flag,est_gamma_parm,patch_list_file,stamps_PART_limitation)
%STAMPS Stanford Method for Persistent Scatterers
%   STAMPS(START_STEP,END_STEP,PATCHES_FLAG,EST_GAMMA_FLAG) Default is to run all steps.
%   A subset of steps may be selected with START_STEP and/or END_STEP
%   STEP 1 = Initial load of data
%   STEP 2 = Estimate gamma 
%   STEP 3 = Select PS pixels
%   STEP 4 = Weed out adjacent pixels
%   STEP 5 = Correct wrapped phase for spatially-uncorrelated look angle error and merge patches
%   STEP 6 = Unwrap phase
%   STEP 7 = Calculate spatially correlated look angle (DEM) error 
%   STEP 8 = Filter spatially correlated noise 
%   STEP 0 = Continue from the last known stage till the end-stage selected
%   
%   PATCHES_FLAG Default 'y'. Set to 'n' to process all data as one patch
%
%   EST_GAMMA_PARM is an optional parameter passed to PS_EST_GAMMA_QUICK
%
%   PATCH_LIST_FILE is an optional argument specifying the file list of
%   patches to be processed. Note that from step 5 and above one should use
%   all patches to merge results.
%
%   If current directory is a single patch, stamps only operates in the
%   current directory, but if current directory contains many patches,
%   stamps operates on them all.
%
%   Andy Hooper, June 2006

%   =================================================================
%   07/2006 AH: END_STEP added
%   09/2006 AH: ps_load removed (obsolete)
%   09/2006 AH: small baselines added 
%   11/2006 AH: patches added
%   01/2007 AH: calculate spatially correlated look angle error added
%   03/2009 AH: simultaneously estimate velocity when SCLA estimated
%   03/2009 AH: smooth SCLA for unwrapping iteration
%   03/2010 AH: move ps_cal_ifg_std to after merge step
%   12/2012 AH: add gamma option
%   12/2012 DB: add patch_list_file argument as option
%   09/2013 DB: update the stamps version number 
%   09/2015 DB: Check if patches do have PS before proceeding with
%               processing.
%   09/2015 DB: Fix when running stamps in a patch folder mode when no PS are left
%   09/2015 AH: allow for non-differentiation of caps by dir
%   01/2016 DB: include stamps_save in step 1-4.
%   08/2016 AH: Fix bug of scn_kriging_flag not being set
%   06/2017 DB: Catching when no PS are left from step 1, allow for re-run
%               when parameters have changed.
%   06/2017 DB: Option to continue from last know processing step
%   08/2107 AH: Removed catch as proceeds also when error in Step 1
%   =================================================================

nfill=40;
fillstr=[repmat('#',1,nfill),'\n'];
skipstr='\n';
msgstr=fillstr;

fprintf(skipstr);
logit(fillstr);
msgstr(round(nfill)/2-12:round(nfill/2)+13)=' StaMPS/MTI Version 4.0b6 ';
logit(msgstr);
msgstr(round(nfill)/2-12:round(nfill/2)+13)='  Beta version, Jun 2018  ';
logit(msgstr);
logit(fillstr);
fprintf(skipstr);



quick_est_gamma_flag=getparm('quick_est_gamma_flag');
reest_gamma_flag=getparm('select_reest_gamma_flag');
unwrap_method=getparm('unwrap_method');
unwrap_prefilter_flag=getparm('unwrap_prefilter_flag');
small_baseline_flag=getparm('small_baseline_flag');
insar_processor=getparm('insar_processor');
scn_kriging_flag=getparm('scn_kriging_flag');

if nargin<1 || isempty(start_step)==1
    start_step=1;
end

if nargin<2 || isempty(end_step)==1
    end_step=8;
end

if nargin<3 || isempty(patches_flag)==1
    if start_step<6
        patches_flag='y';
    else
        patches_flag='n';
    end
end

if nargin<4 || isempty(est_gamma_parm)==1
    est_gamma_parm=0;
end

if nargin<5 || isempty(patch_list_file)     % [DB] allow for own specified patch list file
    patch_list_file = 'patch.list';
    new_patch_file = 0;
else
    % use own file
    new_patch_file = 1;
end

% In support of the multi-core option limit processing to steps 1-5a,
% 5b-upwards, or the old method where all is processed together.
if nargin<6 || isempty(stamps_PART_limitation)
    stamps_PART_limitation=0;
end
stamps_PART1_flag='y';
stamps_PART2_flag='y';
if stamps_PART_limitation==1
    stamps_PART2_flag='n';
end
if stamps_PART_limitation==2
    stamps_PART1_flag='n';
end

if strcmpi(patches_flag,'y')
    if exist(patch_list_file,'file')
        fid=fopen(patch_list_file);
        i=0;
        while 1
            nextline=fgetl(fid);
            if ischar(nextline)
                i=i+1;
                patchdir(i).name=nextline;
            else
                break
            end
        end
	fclose(fid);
    else
        patchdir=dir('PATCH_*');
        patchdir = patchdir(find(~cellfun(@(x) strcmpi(x,'patch_noover.in'),{patchdir(:).name})));
    end
    if isempty(patchdir)
        patches_flag='n';
    else
        ps_parms_default
        patches_flag='y';
    end
end

if ~strcmpi(patches_flag,'y')
    patchdir(1).name='.';
    logit('Will process current directory only')
else
    logit('Will process patch subdirectories')
end

currdir=pwd;

nfill=40;
fillstr=[repmat('#',1,nfill),'\n'];
msgstr=fillstr;





% limit the processing to step 1-5a
start_step_or = start_step;
if strcmpi(stamps_PART1_flag,'y')
  for i=1:length(patchdir)
    if ~isempty(patchdir(i).name)
      cd(patchdir(i).name)
      patchsplit=strsplit(pwd,'/');
      %fprintf(skipstr);
      %logit(sprintf('Processing %s',patchsplit{end}))
    
      % store if patch dir is empty
      if exist('no_ps_info.mat','file')~=2
         stamps_step_no_ps = zeros([5 1 ]);       % keep for the first 5 steps only
         save('no_ps_info.mat','stamps_step_no_ps')
      end
      

        % if start_step is 0, then start from the latest stage it was
        if start_step_or==0
            % check the processing stage of stamps for all the patches
            % step 4 find a ps_weed file
            % step 3 find a ps_select file
            % step 2 find a pm file
            % step 1 find a ps file
            % or no PS in case stamps_step_no_ps is found with a 1 in there 

            if exist('weed1.mat','file')==2
                   start_step=5;
                   setpsver(2);
            elseif exist('select1.mat','file')==2
                   start_step=4;
            elseif exist('pm1.mat','file')==2
                   start_step=3;
            elseif exist('ps1.mat','file')==2
                   start_step=2;
            else
                start_step=1;
            end

            if start_step>end_step
                fprintf(['\n' patchsplit{end} ': already up to end stage ' num2str(end_step) ' \n'])
            else
                fprintf(['\n' patchsplit{end} ': complete up to stage ' num2str(start_step-1) ' \n'])
            end
        end



      
      if start_step==1
        msgstr(round(nfill)/2-3:round(nfill/2)+4)=' Step 1 ';
        fprintf(skipstr);
        logit(fillstr);
        logit(msgstr);
        logit(fillstr)
        logit(['Directory is ',patchsplit{end}])
        fprintf(skipstr);
        if strcmpi(small_baseline_flag,'y')
%             try 
                if strcmpi(insar_processor,'gamma') | strcmpi(insar_processor,'snap')
                    sb_load_initial_gamma;
                elseif strcmpi(insar_processor,'gsar')
                    sb_load_initial_gsar;
                elseif  strcmpi(insar_processor,'isce')
                    if exist('data_inc','var')==0
                        % already in patch dir, file contained in the InSAR dir
                        inc_angle = ['..' filesep 'inc_angle.raw'];
                        if exist(inc_angle,'file')~=2
                             inc_angle = ['..' filesep inc_angle];
                        end
                        if exist(inc_angle,'file')==2
                            fprintf('Found inc angle file, will load the data \n')
                            data_inc = (load_isce(inc_angle));
                        else
                            data_inc=[];
                        end
                    end
                    sb_load_initial_isce(data_inc)
                else
                    sb_load_initial;
                end
                load('no_ps_info.mat');
                % reset as we are currently re-processing
                stamps_step_no_ps(1:end)=0;
                
%             catch
% 
%                load('no_ps_info.mat');
%                % reset as we are currently re-processing
%                stamps_step_no_ps(1:end)=0;
%                fprintf('***No PS points left. Updating the stamps log for this****\n')
%                % update the flag indicating no PS left in step 1
%                stamps_step_no_ps(1)=1;
%                psver =1;
%                save('psver.mat','psver')
%                 
%             end
            save('no_ps_info.mat','stamps_step_no_ps')

        else
%             try 
                if strcmpi(insar_processor,'gamma') | strcmpi(insar_processor,'snap')
                    ps_load_initial_gamma;
                elseif strcmpi(insar_processor,'gsar')
                    ps_load_initial_gsar;
                elseif  strcmpi(insar_processor,'isce')
                     if exist('data_inc','var')==0
                        % already in patch dir, file contained in the InSAR dir
                        inc_angle = ['..' filesep 'inc_angle.raw'];
                        if exist(inc_angle,'file')~=2
                             inc_angle = ['..' filesep inc_angle];
                        end
                        if exist(inc_angle,'file')==2
                            fprintf('Found inc angle file, will load the data \n')
                            data_inc = (load_isce(inc_angle));
                        else
                            data_inc=[];
                        end
                    end
                    ps_load_initial_isce(data_inc)  
         
                else
                    ps_load_initial;
                end
                load('no_ps_info.mat');
                % reset as we are currently re-processing
                stamps_step_no_ps(1:end)=0;
%             catch
%                 load('no_ps_info.mat');
%                 % reset as we are currently re-processing
%                 stamps_step_no_ps(1:end)=0;
%                 fprintf('***No PS points left. Updating the stamps log for this****\n')
%                 % update the flag indicating no PS left in step 1
%                 stamps_step_no_ps(1)=1;
%                 save('no_ps_info.mat','stamps_step_no_ps')
%                 psver =1;
%                 save('psver.mat','psver')
% 
%             end
            save('no_ps_info.mat','stamps_step_no_ps')
        end
        elseif start_step <=4
            setpsver(1)
      end

        
      
      
        if start_step<=2 & end_step >=2 
            msgstr(round(nfill)/2-3:round(nfill/2)+4)=' Step 2 ';
            fprintf(skipstr);
            logit(fillstr);
            logit(msgstr);
            logit(fillstr)
            logit(['Directory is ',patchsplit{end}])
            fprintf(skipstr);

            % check if step 1 had more than 0 PS points
            load('no_ps_info.mat');
            % reset as we are currently re-processing
            stamps_step_no_ps(2:end)=0;
            
            % run step 2 when there are PS left in step 1
            if stamps_step_no_ps(1)==0
                if strcmpi(quick_est_gamma_flag,'y')
                    ps_est_gamma_quick(est_gamma_parm);
                else
                    ps_est_gamma(est_gamma_parm);
                end
            else
                stamps_step_no_ps(2)=1;
                fprintf('No PS left in step 1, so will skip step 2 \n')
            end  
            save('no_ps_info.mat','stamps_step_no_ps')
        end

        if start_step<=3 & end_step >=3 
            msgstr(round(nfill)/2-3:round(nfill/2)+4)=' Step 3 ';
            fprintf(skipstr);
            logit(fillstr);
            logit(msgstr);
            logit(fillstr)
            logit(['Directory is ',patchsplit{end}])
            fprintf(skipstr);

            
            
            % check if step 2 had more than 0 PS points
            load('no_ps_info.mat');
            % reset as we are currently re-processing
            stamps_step_no_ps(3:end)=0;
            
            % run step 3 when there are PS left in step 2
            if stamps_step_no_ps(2)==0
                if strcmpi(quick_est_gamma_flag,'y') & strcmpi(reest_gamma_flag,'y')
                    ps_select;
                else
                    ps_select(1);
                end
            else
                fprintf('No PS left in step 2, so will skip step 3 \n')
                stamps_step_no_ps(3)=1;
            end              
            save('no_ps_info.mat','stamps_step_no_ps')
        end

        if start_step<=4 & end_step >=4 
            msgstr(round(nfill)/2-3:round(nfill/2)+4)=' Step 4 ';
            fprintf(skipstr);
            logit(fillstr);
            logit(msgstr);
            logit(fillstr)
            logit(['Directory is ',patchsplit{end}])
            fprintf(skipstr);

            % check if step 3 had more than 0 PS points
            load('no_ps_info.mat');
            % reset as we are currently re-processing
            stamps_step_no_ps(4:end) =0;       % keep for the first 5 steps only
            

            % run step 4 when there are PS left in step 3
            if stamps_step_no_ps(3)==0
                if strcmpi(small_baseline_flag,'y')
                    ps_weed(0,1);
                else
                    ps_weed;
                end
            else
                fprintf('No PS left in step 3, so will skip step 4 \n')
                stamps_step_no_ps(4)=1;

            end
            save('no_ps_info.mat','stamps_step_no_ps')
        end

        if start_step<=5 & end_step >=5 
            msgstr(round(nfill)/2-3:round(nfill/2)+4)=' Step 5 ';
            fprintf(skipstr);
            logit(fillstr);
            logit(msgstr);
            logit(fillstr)
            logit(['Directory is ',patchsplit{end}])
            fprintf(skipstr);


            % check if step 4 had more than 0 PS points
            load('no_ps_info.mat');
            % reset as we are currently re-processing
            stamps_step_no_ps(5:end) = 0;       % keep for the first 5 steps only
            
            % run step 5 when there are PS left in step 3
            if stamps_step_no_ps(4)==0
                ps_correct_phase;
            else
                fprintf('No PS left in step 4, so will skip step 5 \n')
                stamps_step_no_ps(5)=1;
            end
            save('no_ps_info.mat','stamps_step_no_ps')
        end


        cd(currdir)
      end
    end
end

patchsplit=strsplit(pwd,'/');


% check if one can process second part of step 5b and above
if strcmpi(stamps_PART2_flag,'y')
    %%% Loop throught the patches and update the patch.list and keep only those
    %%% that have PS left.
    if patches_flag=='y'
        % go in reverse order such patches can be dropped when needed

        fid = fopen('patch.list_new','w');
        for i=1:length(patchdir)
            % check the file with the PS information
            filename_PS_check = [patchdir(i).name filesep 'no_ps_info.mat'];

            % assume by default to keep patch for backward compatibility
            keep_patch = 1;
            if exist(filename_PS_check,'file')==2
                load(filename_PS_check)
                if sum(stamps_step_no_ps)>=1
                   keep_patch=0; 
                end
            end

            % update the patch list.
            if keep_patch==1
                fprintf(fid,[patchdir(i).name '\n']);
            end
            if i==length(patchdir)
               fclose(fid) ;
            end
        end

        % update the files such in futhre the new patch list will be used.
        movefile('patch.list','patch.list_old');
        movefile('patch.list_new','patch.list');
    end


    if start_step<=5 & end_step >=5 
        abord_flag=0;
        if patches_flag=='y'
            fprintf(skipstr);
            logit(['Directory is ',patchsplit{end}])
            fprintf(skipstr);
            ps_merge_patches
        else
            % this is processing of an individual patch
            % see if there are any PS left
            if exist('no_ps_info.mat','file')==2
                load('no_ps_info.mat')
                if sum(stamps_step_no_ps)>=1
                   abord_flag=1; 
                end
            end
        end

        % see if step 5 can be ran
        if abord_flag==0 
            ps_calc_ifg_std;
        else
            fprintf('No PS left in step 4, so will skip step 5 \n')
        end
    end


    if start_step<=6 & end_step >=6 
        msgstr(round(nfill)/2-3:round(nfill/2)+4)=' Step 6 ';
        fprintf(skipstr);
        logit(fillstr);
        logit(msgstr);
        logit(fillstr)
        logit(['Directory is ',patchsplit{end}])
        fprintf(skipstr);

        ps_unwrap
        if strcmpi(small_baseline_flag,'y')
            sb_invert_uw
        end
    end

    if start_step<=7 & end_step >=7 
        msgstr(round(nfill)/2-3:round(nfill/2)+4)=' Step 7 ';
        fprintf(skipstr);
        logit(fillstr);
        logit(msgstr);
        logit(fillstr)
        logit(['Directory is ',patchsplit{end}])
        fprintf(skipstr);

        if strcmpi(small_baseline_flag,'y')
            ps_calc_scla(1,1)   % small baselines
            ps_smooth_scla(1)
            ps_calc_scla(0,1) % single master
        else
            ps_calc_scla(0,1)
            ps_smooth_scla
        end
    end

    if start_step<=8 & end_step >=8
        msgstr(round(nfill)/2-3:round(nfill/2)+4)=' Step 8 ';
        fprintf(skipstr);
        logit(fillstr);
        logit(msgstr);
        logit(fillstr)
        logit(['Directory is ',patchsplit{end}])
        fprintf(skipstr);

        if strcmpi(scn_kriging_flag,'y')
            ps_scn_filt_krig
        else
            ps_scn_filt
        end
    end
end

logit(1);
