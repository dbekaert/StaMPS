function [saveampname]=ps_load_mean_amp()
%PS_LOAD_MEAN_AMP Load mean_amp file into a matlab workspace
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   02/2010 AH: Reduce memory needs by reading in chunks
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
    error('Cannot find mean_amp.flt in this directory or parent')
end

amp_high=[];
fid=fopen(ampname,'r');
while ~feof(fid) % first read through
  amp_bit=fread(fid,[width,1000],'float');
  amp_bit(isnan(amp_bit))=1;
  amp_bit=sort(amp_bit(:));
  amp_high=[amp_high;amp_bit(round(length(amp_bit)*0.99):end)];
end
fclose(fid);
amp_high=sort(amp_high);
amp_max=log(amp_high(round(length(amp_high)/2)));
clear amp_high

fid=fopen(ampname,'r');
amp_mean=[];
while ~feof(fid) % second read through
  amp_bit=fread(fid,[width,1000],'float');
  amp_bit=amp_bit';
  amp_bit(amp_bit<1)=1;
  amp_bit=log(amp_bit);
  amp_bit(isnan(amp_bit))=0;
  amp_bit(amp_bit>amp_max)=amp_max;
  amp_bit=uint8(amp_bit*255/amp_max);
  amp_mean=[amp_mean;amp_bit];
end

save(saveampname,'amp_mean');


