function [patchdir_new,patchdir_fail,patchdir_reject]= patch_search(patch_list_file,hgt_threshold,bbox_SNWE)
% function to search for patches in a BBox or those below a given elevation.
% Assumes zero elevation for default and the complete dataset as BBOX
% Bekaert David     -  June 2017


check_bbox = 'n';
if nargin<1 || isempty(patch_list_file)     % [DB] allow for own specified patch list file
    patch_list_file = 'patch.list';
end
if nargin<2
    hgt_threshold =-inf;
end
if nargin>2 & ~isempty(bbox_SNWE)
    check_bbox = 'y';
end
fprintf(['--------------------------------------------------- \n']);
fprintf(['Get all patches: \t\t\t > ' num2str(hgt_threshold) 'm\n']);
fprintf(['Note 0m might not be water level (depends on location)\n']);
if strcmpi(check_bbox,'y')
    fprintf(['Get all patches within bbox (SNWE): \t ' num2str(bbox_SNWE(1)) ' ' num2str(bbox_SNWE(2)) ' ' num2str(bbox_SNWE(3)) ' ' num2str(bbox_SNWE(4))  '\n']);
end
fprintf(['Outputs will be written in: \t\t patch_new.list \n'])
fprintf(['--------------------------------------------------- \n']);

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

counter = 1;
counter1 = 1;
counter2 = 1;
patchdir_new = [];
patchdir_fail = [];
patchdir_reject = [];
h1 = figure;
if strcmpi(check_bbox,'y')
    hold on
    bbox_SNWE_poly = [bbox_SNWE(3) bbox_SNWE(1) ; bbox_SNWE(3) bbox_SNWE(2) ;bbox_SNWE(4) bbox_SNWE(2)  ;bbox_SNWE(4) bbox_SNWE(1) ;bbox_SNWE(3) bbox_SNWE(1)];
    plot(bbox_SNWE_poly(:,1),bbox_SNWE_poly(:,2),'k-','linewidth',2);
    legend_str{1} = 'region';
end
for i=1:length(patchdir)
    if ~isempty(patchdir(i).name)
      cd(patchdir(i).name)
      patchsplit=strsplit(pwd,'/');
      if exist('hgt1.mat','file')==2
          load('hgt1.mat')
          check_hgt = sum(hgt>hgt_threshold)>0;
          
          % getting lon lat info
          ll = load('ps1.mat','lonlat');
          ll_mean = mean(ll.lonlat,1);
          % if requested find patches within a bbox
          in_box = logical(1);
          if strcmpi(check_bbox,'y')
             ix = inpolygon(ll.lonlat(:,1),ll.lonlat(:,2),bbox_SNWE_poly(:,1),bbox_SNWE_poly(:,2)); 
             in_box = sum(ix)>1;
          end
          % check if the heights are met
          if check_hgt
  
      
             try 
                figure(h1)
                hold on
                if in_box
                    plot(ll_mean(1),ll_mean(2),'go','markerfacecolor','g')
                    patchdir_new{counter}=patchsplit{end};
                    counter = counter+1;
             
                else
                    plot(ll_mean(1),ll_mean(2),'go')
                end
             catch
             end

          else
             patchdir_reject{counter1}=patchsplit{end};
             counter1 = counter1+1; 
             try 
                figure(h1)
                hold on
                plot(ll_mean(1),ll_mean(2),'ro')
             catch
             end
          end
                   
      else
         patchdir_fail{counter2}= patchsplit{end} ;
         counter2 = counter2+1;
      end
      cd ..
      
    end
end
try 
    figure(h1)
    box on
    set(gca,'fontsize',15)
catch
end

% write out new patch_new.list
fid = fopen('patch_new.list','w');
for k=1:length(patchdir_new);
   fprintf(fid,[patchdir_new{k} '\n']) ;
end
fclose(fid);

% write out new patch_fail.list
fid = fopen('patch_fail.list','w');
for k=1:length(patchdir_fail);
   fprintf(fid,['rm -rf ' patchdir_fail{k} '\n']) ;
end
fclose(fid);

% write out new patch_reject_region.list
fid = fopen('patch_reject_region.list','w');
for k=1:length(patchdir_reject);
   fprintf(fid,['rm -rf ' patchdir_reject{k} '\n']) ;
end
fclose(fid);

if length(patchdir_fail)>0
    fprintf('Some patches failes to load and are not plotted, check second output argument! \n')
end



