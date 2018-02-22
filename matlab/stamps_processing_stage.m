function []= stamps_processing_stage(patch_list_file,stage_threshold)
% checks the processing stage of stamps for all the patches and optional
% output generation of those still to run up to stage_threshold.
% step 4 find a ps_weed file
% step 3 find a ps_select file
% step 2 find a pm file
% step 1 find a ps file
% or no PS in case stamps_step_no_ps is found with a 1 in there    => gets value 5
%
% INPUTS:
% stage_threshold                   Patches which have not been processed up to stage_threshold level.
%
% OUTPUTS:
% overview figure with processing stages
% optional a new patch_stage.list with those that still need to be run up to the stagethreshold level
%
% Bekaert David         -  June 2017
% mofidications
% DB    02/2018     Adding patch numbers to the plot

if nargin<1 || isempty(patch_list_file)     % [DB] allow for own specified patch list file
    patch_list_file = 'patch.list';
end
if nargin<2 || isempty(stage_threshold)
    stage_threshold=4;
end
generate_output = 'n';
if nargin >1 
    generate_output = 'y';
end


% actual code
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
    patchdir=dir('PATCH_*');
    patchdir = patchdir(find(~cellfun(@(x) strcmpi(x,'patch_noover.in'),{patchdir(:).name})));
end

% tracking the processing stage
stages = zeros([length(patchdir) 1]);
ll_mean = NaN([length(patchdir) 2]);
counter = 1;
patchID = [];
for i=1:length(patchdir)
    if ~isempty(patchdir(i).name)
      cd(patchdir(i).name)
      patchsplit=strsplit(pwd,'/');
      patchdir_keep{counter}=patchsplit{end};

      patchID = [patchID ; str2num(patchdir_keep{counter}(7:end))];
      
      % stage 4
      if exist('weed1.mat','file')==2
           stages(i,1)=4;
      elseif exist('select1.mat','file')==2
           stages(i,1)=3;
      elseif exist('pm1.mat','file')==2
           stages(i,1)=2;
      elseif exist('ps1.mat','file')==2
           stages(i,1)=1;
      end
      % see if the processing was once terminated
      if exist('no_ps_info.mat','file')==2
            load('no_ps_info.mat');
            ix = find(stamps_step_no_ps==1);
            if ~isempty(ix)
               stages(i,1)=5;
            end
      end
      
      
      try
         % getting lon lat info
          ll = load('ps1.mat','lonlat');
          ll_mean(i,:) = mean(ll.lonlat,1);
      catch
          
      end
      counter = counter+1;
      cd ..
    end
end


figure('name',['Processing ' datestr(clock,'yyyymmmdd-HHMM')]);
scatter3(ll_mean(:,1),ll_mean(:,2),stages,15,stages,'filled');
for k=1:length(patchID)
    if floor(k/10)*10 ==k
        text(ll_mean(k,1),ll_mean(k,2),num2str(patchID(k)))
    end
end
view(0,90)
axis equal
axis tight
cc= colorbar;
colormap(lines(6))
caxis([-0.5 5.5])
box on
xlabel('lon')
ylabel('lat')
title(cc,'processing stage')
title([datestr(clock,'yyyy-mmm-dd HH:MM')])
print(gcf,'-dpng',['processing_stage_' datestr(clock,'yyyymmmdd-HHMM') '.png'])


% provide optional output if requested
if strcmpi(generate_output,'y')
    % finding all the patches that are below a stage threshold
    ix = find(stages<stage_threshold & stages>0);
    % write out new patch_stage.list
    if ~isempty(ix)
        fid = fopen('patch_stage.list','w');
        for k=1:length(ix);
            if k==length(ix)
               fprintf(fid,[patchdir_keep{ix(k)} ]) ;
            else
               fprintf(fid,[patchdir_keep{ix(k)} '\n']) ;
            end
        end
        fclose(fid);
    end
end
