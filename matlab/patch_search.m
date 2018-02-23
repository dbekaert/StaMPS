function [patchdir_new,patchdir_fail,patchdir_reject]= patch_search(patch_list_file,hgt_threshold,bbox_SNWE)
% function to search for patches in a BBox or those below a given elevation.
% Assumes zero elevation for default and the complete dataset as BBOX
% Bekaert David     -  June 2017
% modifications:
% DB    02/2018     Allow for use of raw patches too


check_bbox = 'n';
if nargin<1 || isempty(patch_list_file)     % [DB] allow for own specified patch list file
    patch_list_file = 'patch.list';
end
if nargin<2 || isempty(hgt_threshold)
    hgt_threshold =-inf;
end
if nargin>2 & ~isempty(bbox_SNWE)
    check_bbox = 'y';
end

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

% chcking if this needs to be ran on the raw data or now
% check if any stmaps processign has been ran or no
patchdir_hgt = [];
for i=1:length(patchdir)
    if ~isempty(patchdir(i).name)
      cd(patchdir(i).name)
      if exist('hgt1.mat','file')
          patchdir_hgt = [patchdir_hgt ; 1];
      end
      cd ..
    end
end

raw_data_flag = 'n';
if isstr(hgt_threshold)
    raw_data_flag = 'y';
end
if isempty(patchdir_hgt) & strcmpi(raw_data_flag,'n')
    fprintf('Looks like you have not yet ran stamps step 1 will use the raw information\n')
    raw_data_flag = 'y';

end

    
if strcmpi(raw_data_flag,'n')
    fprintf(['--------------------------------------------------- \n']);
    fprintf(['Get all patches: \t\t\t > ' num2str(hgt_threshold) 'm\n']);
    fprintf(['Note 0m might not be water level (depends on location)\n']);
    if strcmpi(check_bbox,'y')
        fprintf(['Get all patches within bbox (SNWE): \t ' num2str(bbox_SNWE(1)) ' ' num2str(bbox_SNWE(2)) ' ' num2str(bbox_SNWE(3)) ' ' num2str(bbox_SNWE(4))  '\n']);
    end
    fprintf(['Outputs will be written in: \t\t patch_new.list \n'])
    fprintf(['--------------------------------------------------- \n']);
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
      
      
      if strcmpi(raw_data_flag,'n')
      
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
          
      else
          % load the lonlat information only once
          if i==1
              lon_file = '../lon.raw';
              lat_file = '../lat.raw';
              if exist(lon_file,'file')~=2
                  lon_file = ['..' filesep lon_file];
              end
              if exist(lat_file,'file')~=2
                  lat_file = ['..' filesep lat_file];
              end
              % load the data
              if exist(lat_file,'file')~=2 & exist(lon_file,'file')~=2
                  error('Could not find the lon lat files')
              end
              lon = single(load_isce(lon_file));
              lat = single(load_isce(lat_file));
              lon(lon==0)=nan;
              lat(lat==0)=nan;
          end
          
          location = load('patch_noover.in');
          lines = [location(1):location(2)];
          pixels = [location(3):location(4)];
          ll.lonlat = [reshape(lon(lines,pixels),[],1) reshape(lat(lines,pixels),[],1) ];
          drop = sum(isnan(ll.lonlat),2)>=1;
          ll.lonlat(drop,:)=[];
          ll_mean = double(mean(ll.lonlat,1));
            
          % if requested find patches within a bbox
          in_box = logical(1);
          if strcmpi(check_bbox,'y')
             ix = inpolygon(ll.lonlat(:,1),ll.lonlat(:,2),bbox_SNWE_poly(:,1),bbox_SNWE_poly(:,2)); 
             in_box = sum(ix)>1;
          end
          
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
             hold on
             temp = patchsplit(end); 
             patchnumber = strtrim(temp{1}(7:end));
             if  floor(str2num(patchnumber)/50)*50==str2num(patchnumber)
                 text(ll_mean(1)+0.001,ll_mean(2)+0.001,patchnumber,'k')
             end
             
          catch
          end
          cd ..
                        
      end
      
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



