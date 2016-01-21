function [saveampname]=ps_load_mean_amp()
%PS_LOAD_MEAN_AMP Load mean_amp file into a matlab workspace
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   02/2010 AH: Reduce memory needs by reading in chunks
%   11/2010 AH: Do amplitude merge here instead of ps_merge_patches.m 
%   11/2010 MA: Fix when length(amp_bit) returns zero at the end
%   ======================================================================

if exist('patch.in','file')
    patch=load('patch.in');
    width=patch(2)-patch(1)+1;
else
    widthname='width.txt';
    if ~exist(widthname,'file')
        widthname= ['../',widthname];
    end
    width=load(widthname);
end

ampname=['./mean_amp.flt'];
saveampname=['./amp_mean.mat'];

if ~exist(ampname,'file')
    ampname=['../mean_amp.flt'];
    saveampname=['../amp_mean.mat'];
end

if ~exist(ampname,'file')
  ampname=['./mean_amp.flt'];
  saveampname=['./amp_mean.mat'];
  fprintf('   Merging mean amplitude files\n')
  widthname='width.txt';
  if ~exist(widthname,'file')
      widthname= ['../',widthname];
  end
  width=load(widthname);
  amp=zeros(width,1,'single');

  if exist('./patch.list','file')
    dirname=struct;
    fid=fopen('patch.list','r');
    i=0;
    while feof(fid)==0
        i=i+1;
        dirname(i).name=fgetl(fid);
    end
    fclose(fid);
  else
    dirname=dir('PATCH_*');
  end
  n_patch=length(dirname);
  for i=1:n_patch
   if ~isempty(dirname(i).name)
    cd(dirname(i).name);
    if exist('./mean_amp.flt','file')    
        pxy=load('patch.in');
        fid=fopen('mean_amp.flt');
        amp(pxy(1):pxy(2),pxy(3):pxy(4))=fread(fid,[pxy(2)-pxy(1)+1,inf],'float=>single');
        fclose(fid);
    end 
    cd ..
   end
  end

  fid=fopen('./mean_amp.flt','w');
  fwrite(fid,amp,'float');
  fclose(fid);
  clear amp
    
%error('Cannot find mean_amp.flt in this directory or parent')
end

amp_high=[];
amp_low=[];
fid=fopen(ampname,'r');
while ~feof(fid) % first read through
  amp_bit=fread(fid,[width,1000],'float');
  amp_bit(isnan(amp_bit))=1;
  amp_bit=sort(amp_bit(:));
  if length(amp_bit)==0  % MA fix
    break                % length(amp_bit) returns zero at the end
  end
  amp_high=[amp_high;amp_bit(round(length(amp_bit)*0.91):end)];
  amp_low=[amp_low;amp_bit(round(length(amp_bit)*0.01):end)];
end
fclose(fid);

amp_low=min(amp_low);
amp_high=sort(amp_high);
amp_min=min(amp_low);
amp_max=log(amp_high(ceil(length(amp_high)/2))-amp_min);
clear amp_high amp_low

fid=fopen(ampname,'r');
amp_mean=[];
while ~feof(fid) % second read through
  amp_bit=fread(fid,[width,1000],'float');
  amp_bit=amp_bit'-amp_min;
  amp_bit(amp_bit<1)=1;
  amp_bit=log(amp_bit);
  amp_bit(isnan(amp_bit))=0;
  amp_bit(amp_bit>amp_max)=amp_max;
  amp_bit=uint8(amp_bit*255/amp_max);
  amp_mean=[amp_mean;amp_bit];
  if length(amp_bit)==0  % MA fix
  break
  end
end

save(saveampname,'amp_mean');


