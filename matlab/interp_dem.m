function []=interp_dem(width,infile,outfile,length,n_interp)
%INTERP_DEM Interpolate DEM in azimuth
%   INTERP_DEM(WIDTH,INFILE,OUTFILE,LENGTH,N_INTERP) WIDTH is width of dem in pixels
%
%   Andy Hooper, Mar 2005
%
% ======================================================
% 08/06 AH: generalized to allow variable filenames
% 11/06 AH: interp used instead of interp2 - less memory
% 05/07 AH: more changes to use less memory
% 10/08 AH: allow for different interpolation factors 
% ======================================================

if nargin < 2
   infile='refdem_1l.raw';
end

if nargin < 3
   outfile='refdem_1li.raw';
end

if nargin < 4
   length=0;
end

if nargin < 5
   n_interp=5;
end

fid=fopen(infile); 
fseek(fid,0,1);
fbytes=ftell(fid);
fseek(fid,0,-1);
if length==0
    pbytes=4
else
    pbytes=fbytes/length/width;
end
if pbytes==4
    fprintf('file looks to be single precision\n')
    dem=fread(fid,[width,inf],'float=>single');
    fclose(fid);
elseif pbytes==8
    fprintf('file looks to be double precision\n')
    dem=fread(fid,[width,inf],'double=>single');
    fclose(fid);
else
    fprintf('file length is %d\n',fbytes)
    fprintf('Number of lines (multilooked) is %d\n',length)
    fprintf('Number of pixels (multilooked) is %d\n',width)
    error('file length is apparently incorrect')
end

dem=dem';
[ny,nx]=size(dem);
ix=find(isnan(dem));
dem_ix=logical(ones(size(dem)));
dem_ix(ix)=0;
demi1=dem;
[X_patch,Y_patch]=meshgrid([1:11],[1:11]);

for i=ix'
    i2=ceil(i/ny);
    i1=i-(i2-1)*ny;
    ix1=max(1,i1-5):min(i1+5,ny);
    ix2=max(1,i2-5):min(i2+5,nx);
    dem_patch=dem(ix1,ix2);
    ix_patch=dem_ix(ix1,ix2);
    if sum(ix_patch(:))<40
        ix_patch(:)=1;
    end
    [X_patch,Y_patch]=meshgrid([1:length(ix2)],[1:length(ix1)]);
    demi1(i)=griddata(X_patch(ix_patch),Y_patch(ix_patch),dem_patch(ix_patch),i2-ix2(1)+1,i1-ix1(1)+1,'linear',{'QJ'});
end
demi1(isnan(demi1))=0;
demi2=zeros(size(demi1,1)*n_interp,size(demi1,2));

for i=1:size(demi1,2)
    demi2(:,i)=interp(demi1(:,i),n_interp); % oversample
end

demi2=single(demi2);
demi2=[demi2(1,:);demi2(1,:);demi2(1:end-2,:)];

fid=fopen(outfile,'w');
fwrite(fid,demi2','float');
fclose(fid);
