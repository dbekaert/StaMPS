function []=writecpx(fname,vname,precision)

if nargin < 2
  error('syntax is writecpx(FILE_NAME, VAR_NAME, PRECISION)')
end

if nargin<3
  precision='float';
end

fid=fopen(fname,'w');

vname_flt=zeros(size(vname,1),size(vname,2)*2);
vname_flt(:,1:2:end)=real(vname);
vname_flt(:,2:2:end)=imag(vname);

fwrite(fid,vname_flt.',precision);



fclose(fid);
