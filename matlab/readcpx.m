function [vname]=readcpx(fname,width,lines,endian,skipbytes,precision)
%READCPX Read complex interleaved files
%   A = READCPX(FILE_NAME, WIDTH, LINES, ENDIAN, SKIPBYTES, PRECISION)
%       FILE_NAME and WIDTH are compulsory
%
%   Andy Hooper, Nov 2003

if nargin < 2
  error('syntax is A = readcpx(FILE_NAME, WIDTH, LINES, ENDIAN, SKIPBYTES, PRECISION)')
end

if nargin < 3
    lines = inf;
end

if nargin < 4
    endian = 'n';
end

if nargin < 5
    skipbytes = 0;
end

if nargin < 6
    precision = 'float';
end

fid=fopen(fname,'r',endian);

fseek(fid,skipbytes,-1);

vname=fread(fid,[width*2,lines],[precision,'=>single']).';

vname=(complex(vname(:,1:2:end),vname(:,2:2:end)));
