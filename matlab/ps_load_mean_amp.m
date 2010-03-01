function [saveampname]=ps_load_mean_amp()
%PS_LOAD_MEAN_AMP Load mean_amp file into a matlab workspace
%
%   Andy Hooper, June 2006

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

fid=fopen(ampname,'r');
amp_mean=fread(fid,[width,inf],'float');
amp_mean=amp_mean';
fclose(fid);

amp_mean(amp_mean<1)=1;
amp_mean=log(amp_mean);
amp_mean(isnan(amp_mean))=0;
aa=sort(amp_mean(:));
i_max=round(length(aa)*0.995);
amp_mean(amp_mean>aa(i_max))=aa(i_max);
amp_mean=uint8(amp_mean*255/max(amp_mean(:)));
save(saveampname,'amp_mean')

