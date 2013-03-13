function stamps_mc_header(start_step,end_step,est_gamma_parm,n_cores,patch_list_file)
% stamps_mc_header(start_step,end_step,est_gamma_parm,n_cores,patch_list_file)
% Program mainly for the first 5 steps of Stamps that splits up the
% parameter list based on the number of processing cores specified. An
% upper bound is set by the cores available.
% INPUTS:
%
%
% By David Bekaert - PhD student - University of Leeds
% December 2012
% modifications:
%


if nargin<3 || isempty(est_gamma_parm)
   est_gamma_parm=[];
end
if nargin<4 || isempty(n_cores)
    n_cores=1;
end
if nargin<5 || isempty(patch_list_file)
   patch_list_file='patch.list';
end
patches_flag=[];

% checking the number of cores by connecting them once all
% matlabpool('open')
% n_cores_available = matlabpool('size');
% matlabpool('close');
n_cores_available = 12;     % maximum for local setup machine

if n_cores_available<n_cores
   fprintf(['Too many cores, maximum cores available: ' num2str(n_cores_available) ' ... \n'])
   n_cores = n_cores_available;
end
if n_cores_available==n_cores
   fprintf('No others cores free once launched ... \n') 
end


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
    n_patches=1;
end
% checking the number of patches to be processed per core
if n_patches<n_cores
    fprintf('Decrease number of cores, as lesser patches are being processed... \n')
   n_cores =  n_patches;
end
n_patches_core = ceil(n_patches./n_cores);
% splitting of the patch list
for k=1:n_cores
    ix = [(k-1)*n_patches_core+1:1:k*n_patches_core];
    ix(ix>n_patches) =[];
    patch_list_filename = ['patch.list_split_' num2str(k)];
    fid = fopen(patch_list_filename,'w');
    for ll=1:length(ix)
        str_temp = patchdir(ix(ll));
        str_temp = str_temp.name;
        copyfile('parms.mat',[str_temp filesep]);
        fprintf(fid,'%s\n',str_temp);
    end
    fclose(fid);
end

keyboard

% openign the set of processors which will be used
matlabpool('open',n_cores)
% launch of the processing jobs
parfor k=1:n_cores
    stamps_mc(start_step,end_step,patches_flag,est_gamma_parm,['patch.list_split_' num2str(k)])
end
% closing the pool of workers again
matlabpool('close')



% Checking if the processing went fine for each core


