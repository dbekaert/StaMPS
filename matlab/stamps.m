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
%
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
%   =================================================================

nfill=40;
fillstr=[repmat('#',1,nfill),'\n'];
skipstr='\n';
msgstr=fillstr;

msgstr(round(nfill)/2-11:round(nfill/2)+12)=' StaMPS/MTI Version 3.3 ';
fprintf([skipstr,fillstr,msgstr]);
msgstr(round(nfill)/2-11:round(nfill/2)+12)=' Beta version, Sep 2013 ';
fprintf([msgstr,fillstr,skipstr]);


quick_est_gamma_flag=getparm('quick_est_gamma_flag');
unwrap_method=getparm('unwrap_method');
unwrap_prefilter_flag=getparm('unwrap_prefilter_flag');
small_baseline_flag=getparm('small_baseline_flag');
insar_processor=getparm('insar_processor');


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
    else
        patchdir=dir('PATCH*');
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
    fprintf('Processing current directory only...\n')
else
    fprintf('Processing patch directories...\n')
end

currdir=pwd;

nfill=40;
fillstr=[repmat('#',1,nfill),'\n'];
msgstr=fillstr;


% limit the processing to step 1-5a
if strcmpi(stamps_PART1_flag,'y')
    for i=1:length(patchdir)
        cd(patchdir(i).name)
        logit(sprintf('\nProcessing PATCH %s',pwd))

        if start_step==1
            msgstr(round(nfill)/2-7:round(nfill/2)+7)=' StaMPS Step 1 ';
            fprintf([skipstr,fillstr,msgstr,fillstr,skipstr]);
            if strcmpi(small_baseline_flag,'y')
                if strcmpi(insar_processor,'gamma')
                    sb_load_initial_gamma;
                else
                    sb_load_initial;
                end
            else
                if strcmpi(insar_processor,'gamma')
                    ps_load_initial_gamma;
                else
                    ps_load_initial;
                end
            end
        elseif start_step <=4
            setpsver(1)
        end

        if start_step<=2 & end_step >=2 
            msgstr(round(nfill)/2-7:round(nfill/2)+7)=' StaMPS Step 2 ';
            fprintf([skipstr,fillstr,msgstr,fillstr,skipstr]);
            if strcmpi(quick_est_gamma_flag,'y')
                ps_est_gamma_quick(est_gamma_parm);
            else
                ps_est_gamma(est_gamma_parm);
            end
        end

        if start_step<=3 & end_step >=3 
            msgstr(round(nfill)/2-7:round(nfill/2)+7)=' StaMPS Step 3 ';
            fprintf([skipstr,fillstr,msgstr,fillstr,skipstr]);
            if strcmpi(quick_est_gamma_flag,'y')
                ps_select;
            else
                ps_select(1);
            end
        end

        if start_step<=4 & end_step >=4 
            msgstr(round(nfill)/2-7:round(nfill/2)+7)=' StaMPS Step 4 ';
            fprintf([skipstr,fillstr,msgstr,fillstr,skipstr]);

            % check if step 3 had more than 0 PS points
            if exist('no_ps_info.mat','file')==2
               load('no_ps_info.mat');
               % reset as we are currently re-processing
            else
               % keep backward compatibility
               stamps_step_no_ps = zeros([5 1 ]);       % keep for the first 5 steps only
            end


            % run step 4 when there are PS left in step 3
            if stamps_step_no_ps(3)==0
                if strcmpi(small_baseline_flag,'y')
                    ps_weed(0,1);
                else
                    ps_weed;
                end
            else
                fprintf('No PS left in step 3, so will skip step 4 \n')
            end
        end

        if start_step<=5 & end_step >=5 
            msgstr(round(nfill)/2-7:round(nfill/2)+7)=' StaMPS Step 5 ';
            fprintf([skipstr,fillstr,msgstr,fillstr,skipstr]);


            % check if step 4 had more than 0 PS points
            if exist('no_ps_info.mat','file')==2
               load('no_ps_info.mat');
               % reset as we are currently re-processing
            else
               % keep backward compatibility
               stamps_step_no_ps = zeros([5 1 ]);       % keep for the first 5 steps only
            end

            % run step 5 when there are PS left in step 3
            if stamps_step_no_ps(4)==0
                ps_correct_phase;
            else
                fprintf('No PS left in step 4, so will skip step 5 \n')
            end
        end


        cd(currdir)
    end
end


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
        msgstr(round(nfill)/2-7:round(nfill/2)+7)=' StaMPS Step 6 ';
        fprintf([skipstr,fillstr,msgstr,fillstr,skipstr]);
        ps_unwrap
        if strcmpi(small_baseline_flag,'y')
            sb_invert_uw
        end
    end

    if start_step<=7 & end_step >=7 
        msgstr(round(nfill)/2-7:round(nfill/2)+7)=' StaMPS Step 7 ';
        fprintf([skipstr,fillstr,msgstr,fillstr,skipstr]);
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
        msgstr(round(nfill)/2-7:round(nfill/2)+7)=' StaMPS Step 8 ';
        fprintf([skipstr,fillstr,msgstr,fillstr,skipstr]);
        if strcmpi(scn_kriging_flag,'y')
            ps_scn_filt_krig
        else
            ps_scn_filt
        end
    end
end

