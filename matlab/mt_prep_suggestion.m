function [] = mt_prep_suggestion(small_baseline_flag)
% suggest the user on the amount of patcehs to use
% give as input if you want the suggestion for the SM of SB processing.
% By Bekaert David - Jet Propulsion Laboratory
% Jan 2016
% modifications:
%


if nargin<1
    error('Indicate at least the small_baseline_flag, ''n'' for PS and ''y'' for SB')
end

% see if we are in the SB processing directory
if strcmpi(small_baseline_flag,'y')
    if exist([pwd filesep 'SMALL_BASELINES'],'dir')==7
        cd('SMALL_BASELINES')
    end
end

% get the number of interferograms
if strcmpi(small_baseline_flag,'y')
    fprintf('Getting number of interferograms from SB processing path:\n')
    [a,b] = system(['\ls -d [1,2]* | sed ''' 's/_/ /''' ' > small_baselines_considered.list']);
    temp = load('small_baselines_considered.list');
    temp = [num2str(temp(:,1)) repmat('_',size(temp,1),1) num2str(temp(:,2))];
    n_ifg = size(temp,1);
else
    fprintf('Getting number of interferograms from SM processing path:\n')
    [a,b] = system(['\ls -d [1,2]*  > singlemaster_baselines_considered.list']);
    temp = num2str(load('singlemaster_baselines_considered.list'));
    n_ifg = size(temp,1);
end
fprintf(['number of ifgs = ' num2str(n_ifg) '\n']);

% get the size of the interferograms
ifg_path = [temp(1,:) filesep 'interferogram.out'];
clear temp
fprintf(['Getting interferogram size from : ' ifg_path '\n'])
[a,b]= system(['echo `grep ''Number\ of\ lines\ (multilooked):'' ' ifg_path ' | awk ''NR==1{print $5}''` > templines.txt' ]);
n_lines = load('templines.txt');
fprintf(['number of lines = ' num2str(n_lines) '\n'])
[a,b]= system(['echo `grep ''Number\ of\ pixels\ (multilooked):'' ' ifg_path ' | awk ''NR==1{print $5}''`  > temppixels.txt']);
n_pixels = load('temppixels.txt');
fprintf(['number of pixels = ' num2str(n_lines) '\n'])


% stamps recomendation is to have less than 5 million pixels per patch per SLC.
% with the larger amoutn of itnerferograms this actually becomes tricky to
% define for an SLC. So lets assume that that would have been for a 
% set of 50 interferograms.
n_pixels_suggested_StaMPS = 5*10^6;
n_ifgs_suggested_StaMPS = 50;
n_total_suggested_StaMPS=  n_pixels_suggested_StaMPS*n_ifgs_suggested_StaMPS;

% number of actual pixels in pur dataset
n_total = n_pixels*n_lines*n_ifg;
% number of patches
n_total_patches = ceil(n_total./n_total_suggested_StaMPS);
% lets now define in such a way they are equally split in azimuth and rangedirection


lines_ratio = n_lines./n_pixels;
n_patches_range = round(sqrt(n_total_patches/lines_ratio));
n_patches_azimuth = round(n_patches_range*lines_ratio);
if (n_patches_azimuth==0)
    n_patches_azimuth=1;
end
if (n_patches_range==0)
    n_patches_range=1;
end
fprintf(['\nSuggested mtprep call:\n'])
if strcmpi(small_baseline_flag,'y')
    fprintf(['mt_prep 0.6 ' num2str(n_patches_range) ' ' num2str(n_patches_azimuth) ' 50 200 \n'])
else
    fprintf(['mt_prep 0.4 ' num2str(n_patches_range) ' ' num2str(n_patches_azimuth) ' 50 200 \n'])
end
