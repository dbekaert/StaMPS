function [] = sb_fix_unwrap_manual(ifg_list,reset_flag)
% [] = sb_fix_unwrap_manual(reset_flag)
% script to manually fix unwrapping errors in Small Baseline processing
%
% This function allows to manually fix unwrapping errors by iterating
% the results of ps_plot('rsb'). You can add add or subtract an integer number
% of 2pi cycles over a selected region. You can keep iterating till the 'rsb'
% errors have dissapeared. If it goes wrong and you do want to reset your 
% unwrapped data back to that of the original you can set the reset_flag to 1;
%
% By Bekaert David - University of Leeds - Sept 2014
% 
% modifications:
%   Bekaert David   05/2015     Fix in command line output, fix the code
%                               for ps_plot
% Bekaert David     09/2015     Add file restoring.
% Bekaert David     12/2015     Also plot the unwrapped phase
%

% flags
if nargin<1
    ifg_list =[];
end
if nargin<2
    reset_flag=0;
end
deramp_flag =0;         % optional deramp rsb. Might be needed when processed with gamma

% keeping the original data, incase user want to reset it
if exist('phuw_sb2_original.mat','file')~=2
    copyfile('phuw_sb2.mat','phuw_sb2_original.mat')
    copyfile('phuw2.mat','phuw2_original.mat')
    copyfile('phuw_sb_res2.mat','phuw_sb_res2_original.mat')
else
    copyfile('phuw_sb2.mat','phuw_sb2_previous_temp.mat')
    copyfile('phuw2.mat','phuw2_previous_temp.mat')
    copyfile('phuw_sb_res2.mat','phuw_sb_res2_previous_temp.mat')
end

% use previous or start from scratch
if reset_flag==1
    fprintf('Using the original data\n')
   ph_input = load('phuw_sb2_original');
else
   fprintf('Use data from a previous run\n')
   ph_input = load('phuw_sb2.mat');
end

% loading the interferogram information
ps = load('ps2.mat');

% plotting the current rsb data
if deramp_flag==1
    ps_plot('usb-do',1,0,0,ifg_list);
    set(gcf,'position',[ 75   522   560   420])
    ps_plot('rsb-o',1,0,0,ifg_list);
     set(gcf,'position',[ 1136         515         560         420])
else
    ps_plot('usb-d',1,0,0,ifg_list);
        set(gcf,'position',[ 75   522   560   420])

    ps_plot('rsb',1,0,0,ifg_list);
     set(gcf,'position',[ 1136         515         560         420])
end

% get the interferogram that the user needs to adapt
repeat=1;
while repeat==1
    ix_ifg = input('Which interferogram to you want to correct? ','s');
    ix_ifg = str2num(ix_ifg);
    if isempty(ix_ifg) 
        repeat=1;
    elseif ix_ifg<=ps.n_ifg
        repeat=0;
    else
        fprintf(['Not that many interferograms \n'])
    end
end


% plot the rsb value for this interferogram
% option one can use deramped rsb, this is for teh case teh interferograms
% were not created from relative differences, i.e. each interferogram had a
% baseline estimated and there might be some ramping errors because of that.
if deramp_flag==1
%     ps_plot('usb-do',1,0,0,ix_ifg);
%         set(gcf,'position',[ 75   522   560   420])

    ps_plot('rsb-o',1,0,0,ix_ifg);
     set(gcf,'position',[ 1136         515         560         420])
else
%     ps_plot('usb-d',1,0,0,ix_ifg);
%         set(gcf,'position',[ 75   522   560   420])

    ps_plot('rsb',1,0,0,ix_ifg);
     set(gcf,'position',[ 1136         515         560         420])
end

% gettign the polygon for which the correction will be made
fprintf('Define a polygon by clicking on the figure. Press enter when done. \n')
polygon=ginput;

% loop untill the user is happy with it
continue_flag=1;
while continue_flag==1
    repeat=1;
    while repeat==1
        ix_shift= input('By how many cycles to you want to shift this region? [+-integer]','s');
        ix_shift = str2num(ix_shift);
        if isempty(ix_shift) 
            repeat=1;
        elseif (ix_shift./(round(ix_shift)))~=1
            fprintf(['Needs to be an integer number... \n'])
            repeat =1;
        else
            repeat=0;
        end
    end

    % the shift in radians
    radian_shift = ix_shift*2*pi;


    % modifying the itnerferogram
    ix = inpolygon(ps.lonlat(:,1),ps.lonlat(:,2),polygon(:,1),polygon(:,2));
    ph_uw= ph_input.ph_uw;
    ph_uw(ix,ix_ifg)=ph_uw(ix,ix_ifg)+radian_shift;
    msd = ph_input.msd;

    % saving the data
    save('phuw_sb2.mat','ph_uw','msd');

    % re-running the 
    sb_invert_uw

    if deramp_flag==1
        ps_plot('usb-do',1,0,0,ix_ifg);
            set(gcf,'position',[ 75   522   560   420])

        ps_plot('rsb-o',1,0,0,ix_ifg);
         set(gcf,'position',[ 1136         515         560         420])
    else
        ps_plot('usb-d',1,0,0,ix_ifg);
            set(gcf,'position',[ 75   522   560   420])

        ps_plot('rsb',1,0,0,ix_ifg);
         set(gcf,'position',[ 1136         515         560         420])
    end    
    
    repeat=1;
    while repeat==1
        string= input('retry? [y/n] ','s');
        if strcmpi(string,'y')
            repeat=0;
            continue_flag = 1;
        elseif strcmpi(string,'n')
            repeat=0;
            continue_flag = 0;
            % see if the result needs to be kept or reverted
            repeat2=1;
            while repeat2==1
                action_flag= input('Keep this result [y/n]? ','s');
                if strcmpi(action_flag,'y')
                    repeat2=0;
                elseif strcmpi(action_flag,'n')
                    repeat2=0;
                    
                    % restore the codes
                    copyfile('phuw_sb2_original.mat','phuw_sb2.mat')
                    copyfile('phuw2_original.mat','phuw2.mat')
                    copyfile('phuw_sb_res2_original.mat','phuw_sb_res2.mat')

                else
                    fprintf('y or n ...\n')
                end
            end
        else
            fprintf('y or n ...\n')
        end
    end
    
end
