function stamps_mc_header(start_step,end_step,patches_flag,est_gamma_parm,patch_list_file)
% stamps_mc_header(start_step,end_step,est_gamma_parm,patch_list_file)
% Program can be used for step 1-4s of stamps, and which splits up the
% parameter list based on the number of processing cores specified. 
% incase you have more cores than patches, then the number of cores is
% reduced to the number of patches. The number of cors can be specified by 
% using setparm('n_cores',XX). It is recommended to investigate the logs
% after each step.
%
% INPUTS:
% start_step        Start step, defined the same as for regular stamps.m script
%                   By default this is step 1.
% end_step          End step, defined the same as for regular stamps.m script
%                   By default this is step 4
% patches_flag      Process in patches, by default this is 'y', or leave blank
% est_gamma_parm    est_gamma_parm, by default this is 0, or leave blank
%
%
% GENERATED FILES:
% log_stamps_overview   This file contains the process ID for each launched
%                       matlab job. You can use this to terminate one of
%                       the processed using e.g. top command in linux. The
%                       processes are ranked starting with CORE1 file.
% patch_list_split_XX   The patchlist used in stamps for the XX core.
% log_stamps_split_XX   The matlab command line output for the XX core.
% 
%
% By David Bekaert - PhD student - University of Leeds
% December 2012
% modifications:
% 01/2014   DB  Remove the matlab multi-core functions and use the schell
%               instead by launching multiple matlab jobs
% 01/2014   DB  Integrate steps 1-4 to split patch lists and steps from 5
%               on to use a merged dataset
% 09/2015   DB  Include step 5 in the processing and add some future ideas
% 09/2015   DB  Bug fix 
% 01/2016   DB  Pause between launching cores

% The definition of the stamps steps
if nargin<1 || isempty(start_step)==1
    start_step=1;
end
if nargin<2 || isempty(end_step)==1
    end_step=4;
end
% if end_step>=5
%    fprintf('Multi-core currently only supported step 1-4\n')
%    if start_step>=5
%        error('Will aboard ..., Please proceed with regular stamps.m function \n')
%    end
% end

% keep below such in future there can be a catch to see when step 5 is
% finnished by all cores before proceeding to the second part of step 5.
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
n_cores=getparm('n_cores');


% below converting the input arguments to same strings that can be called
% outside matlab for the multi-core option
if strcmp(patches_flag,'y')
    patches_flag_str='[]';
end
if est_gamma_parm==0
   est_gamma_parm_str='[]';
end
if isnumeric(est_gamma_parm)
    est_gamma_parm_str = num2str(est_gamma_parm);
end

% defining the step_range
step_range = [start_step:end_step];
ix_split_patches = find(step_range<=5);
ix_merged = find(step_range>=5);




% generating a patchlist in case none is given
if ~isempty(ix_split_patches)
    % generating new patch list based on the number of cores selected
    % getting the patch list from the source file
    if exist(patch_list_file,'file')
        fid=fopen(patch_list_file);
        n_patches=0;
        while 1
            nextline=fgetl(fid);
            if ischar(nextline)
                n_patches=n_patches+1;
                patchdir(n_patches).name=nextline;
            else
                break
            end
        end     
    else
        patchdir=dir('PATCH*');
        n_patches=size(patchdir,1);
    end

    
    % checking the number of patches to be processed per core
    if n_patches<n_cores
        fprintf('Decrease number of cores, as lesser patches are being processed... \n')
       n_cores =  n_patches;
    end
    n_patches_core = ceil(n_patches./n_cores);
    % splitting of the patch list
    counter_cores = 0;
    for k=1:n_cores
        ix = [(k-1)*n_patches_core+1:1:k*n_patches_core];
        ix(ix>n_patches) =[];
        patch_list_filename = ['patch_list_split_' num2str(k)];
        for ll=1:length(ix)
            if ll==1
                    fid = fopen(patch_list_filename,'w');
                    counter_cores = counter_cores+1;
            end
            str_temp = patchdir(ix(ll));
            str_temp = str_temp.name;
            if ~isempty(str_temp)
                copyfile('parms.mat',[str_temp filesep]);
                fprintf(fid,'%s\n',str_temp);
            end
            if ll==length(ix)
                fclose(fid);
            end
        end
    end
    n_cores = counter_cores;

    comandstr = 'echo logfile > log_stamps_overview';
    system(comandstr);

    for k=1:n_cores
        comandstr = (['matlab -nodesktop -nodisplay -r "stamps(' num2str(step_range(ix_split_patches(1)) ) ',' num2str(step_range(ix_split_patches(end))) ',' patches_flag_str ',' est_gamma_parm_str ',' '''patch_list_split_' num2str(k) ''',1  ); exit" > log_stamps_split_' num2str(k) ' & ' ]);
        comandstr2 = 'echo "${!}" >> log_stamps_overview';
        command = [comandstr comandstr2];
        [a,b] = system(command);
        fprintf([num2str(k) 'Done, ... pausing 20s for next job launch \n'])
        pause(20)
    end
    
end


% waiting the processors to complete


% run the second component on the merged patches
if ~isempty(ix_merged)
    fprintf(['Once all the processing has completed, i.e. each core, \nrun the following command to complete your processing\n'])
    fprintf(['stamps(' num2str(step_range(ix_merged(1))) ',' num2str(step_range(ix_merged(end))) ',[],0,[],2);\n']);


    % the following can be used in future, but needs a check in the
    % code which waits to proceed untill all cores have finnised
    % processing. See comment above ~L158.
%           stamps(ix_merged(1),ix_merged(end),[],0,[],2);
end




