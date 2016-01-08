function doris_mc_header(operation,n_cores,process_list_file, extra_arguments,extra_arguments_before)
% Program for the StaMPS preprocessing (selected ROI_PAC and DORIS functions)
% doris_mc_header(operation,n_cores,process_list_file,extra_arguments)
% An upper bound is set by the cores available.
% Valid operations are: make_cohs, make_small_baselines,
% make_slcs_envi_vor, make_ifgs, make_resample, make_dems, make_multilook,
% make_coreg_simple, make_coarse, make_slcs_alos
% By David Bekaert - PhD student - University of Leeds
% December 2012
% modifications: 
% 12/12	 DB:	Add multi-looking processing
% 12/12  DB:	Add coherence processing
% 12/12  DB:    Dump output for each core into a screen
% 09/14  DB:    remove matlab multi-core option
% 05/15  TI:    Add 'on the fly' list file construction
% 05/15  DB:    Test and add exceptions
% 12/15  DB:    Adding dem_assist support to multi-core
% 12/15  DB:    Include a pause as last core does not always start
%
% doris_mc_header('make_coarse',20)

% Checking the input arguments
% only allow for the follwoing processing commands
if ~strcmp(operation,'make_cohs') && ~strcmp(operation,'make_small_baselines') && ~strcmp(operation,'make_slcs_envi_vor') && ~strcmp(operation,'make_ifgs') && ~strcmp(operation,'make_resample')  && ~strcmp(operation,'make_dems') && ~strcmp(operation,'make_multilook') && ~strcmp(operation,'make_coreg_simple') && ~strcmp(operation,'make_coarse') && ~strcmp(operation,'make_dem_assist') && ~strcmp(operation,'make_slcs_alos')
    fprintf('not a valid operation to perform \n')
    keyboard
end
% the number of cores
if nargin<2 || isempty(n_cores)
    n_cores=30;
end
disp(['Number of cores used in processing: ',num2str(n_cores)])


if nargin <3
  process_list_file=[];
  disp('Process list file was not given by user.')
end

%% Attempt to generate the input lists on the fly if not specified by the user (files end in tmp)

% N.B. All these operations assume the user is in the correct directory for the process they're attempting to run

if isempty(process_list_file)
    disp('Attempting to construct process list file automatically...')
	% make_cohs (in INSAR_XXXXXXXX folder)
	if strcmp(operation,'make_cohs')
	   system_command = ['ls -d ' pwd '/[1,2]*/interferogram.out | gawk ''BEGIN {FS="interferogram.out"} {print $1}''  > make_cohs.list_tmp'];	
	   [a,b] = system(system_command);
	   process_list_file = 'make_cohs.list_tmp';
    end

	% make_small_baselines (in INSAR_XXXXXXXX folder)
	if strcmp(operation,'make_small_baselines')
	   process_list_file = 'small_baselines.list';
	end

	% make_slcs_envi_vor (in SLC folder)
	if strcmp(operation,'make_slcs_envi_vor')
	   process_list_file = 'make_slcs.list';
	end

	% make_ifgs (in INSAR_XXXXXXXX folder)
	if strcmp(operation,'make_ifgs')
	   system_command = ['ls -d ' pwd '/[1,2]*/coreg.out | gawk ''BEGIN {FS="coreg.out"} {print $1}''  > make_ifgs.list_tmp'];
	   [a,b] = system(system_command);
	   process_list_file = 'make_ifgs.list_tmp';
    end
    
    % make_resample (in INSAR_XXXXXXXX folder)
	if strcmp(operation,'make_resample')
	   system_command = ['ls -d ' pwd '/[1,2]*/coreg.out | gawk ''BEGIN {FS="coreg.out"} {print $1}''  > make_ifgs.list_tmp'];
	   [a,b] = system(system_command);
	   process_list_file = 'make_ifgs.list_tmp';
    end
    
    % make_dems option (in INSAR_XXXXXXXX folder)
	if strcmp(operation,'make_dems')
	   system_command = ['ls -d ' pwd '/[1,2]*/coreg.out | gawk ''BEGIN {FS="coreg.out"} {print $1}''  > make_ifgs.list_tmp'];
	   [a,b] = system(system_command);
	   process_list_file = 'make_ifgs.list_tmp';
    end
    
    % make_multilook (requires extra arguments for width of file (required) and number of range looks (default 4), the aspect ratio (default 5) and the format (default r4))
	% User specified file list is given in the function
	if strcmp(operation,'make_multilook')
        disp('make_multilook requires a user generated list and extra_arguments (but no extra_arguments_before)')
        disp('extra_arguments should consist of at least the width of the file and optionally the number of range looks (assumes 4), aspect ratio (assumes 5) and format (assumes r4)')
	    process_list_file = 'make_multilook.list';
    end
    
    % make_coreg_simple option (in INSAR_XXXXXXXX folder)
	if strcmp(operation,'make_coreg_simple')
	   system_command = ['ls -d ' pwd '/[1,2]*/coreg.out | gawk ''BEGIN {FS="coreg.out"} {print $1}''  > make_ifgs.list_tmp'];
	   [a,b] = system(system_command);
	   process_list_file = 'make_ifgs.list_tmp';
    end

	% make_coarse (in INSAR_XXXXXXXX folder)
	if strcmp(operation,'make_coarse')
	   system_command = ['ls -d ' pwd '/[1,2]*/slave.res | gawk ''BEGIN {FS="slave.res"} {print $1}''  > make_coarse.list_tmp'];
	   [a,b] = system(system_command);
	   process_list_file = 'make_coarse.list_tmp';
    end
    
	% make_dem_assist (in INSAR_XXXXXXXX folder)    
    if strcmp(operation,'make_dem_assist')
	   system_command = ['ls -d ' pwd '/[1,2]*/coreg.out | gawk ''BEGIN {FS="coreg.out"} {print $1}''  > make_dem_assist.list_tmp'];
	   [a,b] = system(system_command);
	   process_list_file = 'make_dem_assist.list_tmp';
	end
    
    % make_slcs_envi_vor (in SLC folder)
	if strcmp(operation,'make_slcs_alos')
	   process_list_file = 'make_slcs.list';
	end

end

% the file containing the paths to be processed
if isempty(process_list_file)
   fprintf('Automatic list file generation failed. You need to specify a list to be processed.');
else
    disp(['List file automatically generated. Using ',process_list_file])
end


if nargin<4
   extra_arguments = ''; 
end
if nargin<5
    extra_arguments_before = '';
end

% generating new file lists based on the number of cores selected
% getting the patch list from the source file
if exist(process_list_file,'file')
    fid=fopen(process_list_file);
    n_jobs=0;
    while 1
        nextline=fgetl(fid);
        if ischar(nextline)
            n_jobs=n_jobs+1;
            patchdir(n_jobs).name=nextline;
        else
            break
        end
    end     
end
% checking the number of patches to be processed per core
if n_jobs<n_cores
    fprintf('Decrease number of cores, as lesser processes are ran than number of cores allocated... \n')
    n_cores =  n_jobs;
end

n_jobs_core = ceil(n_jobs./n_cores);
% splitting of the patch list
for k=1:n_cores
    ix = [(k-1)*n_jobs_core+1:1:k*n_jobs_core];
    ix(ix>n_jobs) =[];
    process_list_filename = [operation '_' process_list_file '_split_' num2str(k)];
    fid = fopen(process_list_filename,'w');
    for ll=1:length(ix)
        str_temp = patchdir(ix(ll));
        str_temp = str_temp.name;
        fprintf(fid,'%s\n',str_temp);
    end
    fclose(fid);
end

fprintf(['The process will be called like: \n' operation ' ' extra_arguments_before ' ' operation '_' process_list_file '_split_' num2str(1) ' ' extra_arguments ' > ' 'CORE_' num2str(1) '_' operation '.log &' '\n'])

fprintf(['If this is correct type the word "return" on the next line and press enter \n'])
keyboard
% launch of the processing jobs
if exist('processor_overview_id.log','file')
    delete(['processor_overview_id.log'])
end
[a,b] = system(['echo files > processor_overview_id.log']);
for k=1:n_cores
    if exist(['CORE_' num2str(k) '_' operation '.log'],'file')
        delete(['CORE_' num2str(k) '_' operation '.log'])
    end
    [a,b] = system([operation ' ' extra_arguments_before ' ' operation '_' process_list_file '_split_' num2str(k) ' ' extra_arguments ' > ' 'CORE_' num2str(k) '_' operation '.log &']);
    [a,b] = system(['echo "${!}" >> processor_overview_id.log']);
    fprintf([num2str(k) ' cores launched \n'])   
    
end

% pausing the system before removing the temp files as some processors
% could be slower
fprintf('Allow all processors to start.\nWill pause the system 60sec before deleting files... \n')
pause(60)

% closing the pool of workers again
for k=1:n_cores
  delete([operation '_' process_list_file '_split_' num2str(k)]);
end

% Remove temporary list files
[a,b] = system(['rm -f *list_tmp']);

