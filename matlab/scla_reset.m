function []=scla_reset(patches_flag)
% SCLA_RESET resets the estimated values of spatially-correlated look angle
%    and master atmosphere/orbit error.
%
%   Andy Hooper, Jan 2007
%
%   ==========================================================
%   03/2009 AH: delete scla_smooth mat files
%   ==========================================================

logit;

if nargin<1
    patches_flag='y';
end

logit(sprintf('Resetting SCLA error and master atmosphere/orbit error...\n'),2)

i=0;

if strcmpi(patches_flag,'y')
    if exist('patch.list','file')
        fid=fopen('patch.list');
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
    end
    if isempty(patchdir)
        patches_flag='n';
    end
end

patchdir(i+1).name='.';


currdir=pwd;

for i=1:length(patchdir)
  if ~isempty(patchdir(i).name)
    cd(patchdir(i).name)
    logit([pwd,' reset'])

    if exist('./psver.mat','file')
      load psver
      sclaname=['scla',num2str(psver),'.mat'];
      if exist([sclaname],'file')
        delete(sclaname)
      end
      sclaname=['scla_smooth',num2str(psver),'.mat'];
      if exist([sclaname],'file')
        delete(sclaname)
      end
      sclaname=['scla_sb',num2str(psver),'.mat'];
      if exist([sclaname],'file')
        delete(sclaname)
      end
      sclaname=['scla_smooth_sb',num2str(psver),'.mat'];
      if exist([sclaname],'file')
        delete(sclaname)
      end
    end

    cd(currdir)
  end
end
logit(1);

