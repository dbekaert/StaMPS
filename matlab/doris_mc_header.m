function doris_mc_header(operation,n_cores,process_list_file, extra_arguments,extra_arguments_before)
% Program for the doris processing
% doris_mc_header(operation,n_cores,process_list_file,extra_arguments)
% An upper bound is set by the cores available.
% Valid operations are: make_cohs make_slcs_envi_vor make_multilook
% By David Bekaert - PhD student - University of Leeds
% December 2012
% modifications: 
% 12/12	 DB:	Add multi-looking processing
% 12/12  DB:	Add coherence processing
% 12/12  DB:    Dump output for each core into a screen
% 09/14  DB:    remove matlab multi-core option



% Checking the input arguments
% only allow for the follwoing processing commands
if ~strcmp(operation,'make_cohs') && ~strcmp(operation,'make_small_baselines') && ~strcmp(operation,'make_slcs_envi_vor') && ~strcmp(operation,'make_ifgs') && ~strcmp(operation,'make_resample')  && ~strcmp(operation,'make_dems') && ~strcmp(operation,'make_multilook') && ~strcmp(operation,'make_coreg_simple') && ~strcmp(operation,'make_coarse') && ~strcmp(operation,'make_slcs_alos')
    fprintf('not a valid operation to perform \n')
    keyboard
end
% the number of cores
if nargin<2 || isempty(n_cores)
    n_cores=30;
end
% the file containign the paths to be processed
if nargin<3 || isempty(process_list_file)
   fprintf('You need to specify a list to be processed');
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

fprintf(['The process will be calles like: \n' operation ' ' extra_arguments_before ' ' operation '_' process_list_file '_split_' num2str(1) ' ' extra_arguments ' > ' 'CORE_' num2str(1) '_' operation '.log &' '\n'])

fprintf(['If this is correct type the word "return" on the next line and press enter \n'])
keyboard
% launch of the processing jobs
delete(['processor_overview_id.log'])
[a,b] = system(['echo files > processor_overview_id.log']);
for k=1:n_cores
    delete(['CORE_' num2str(k) '_' operation '.log'])
    [a,b] = system([operation ' ' extra_arguments_before ' ' operation '_' process_list_file '_split_' num2str(k) ' ' extra_arguments ' > ' 'CORE_' num2str(k) '_' operation '.log &']);
    [a,b] = system(['echo "${!}" >> processor_overview_id.log']);
    fprintf([num2str(k) ' cores launched \n'])   
    
end

% closing the pool of workers again
for k=1:n_cores
 delete([operation '_' process_list_file '_split_' num2str(k)]);
end

