function [] = sb_fix_unwrap_manual(reset_flag,ix_reduced_list)
% [] = sb_fix_unwrap_manual(reset_flag)
% script to manually fix unwrapping errors in Small Baseline processing
%
% This function allows to manually fix unwrapping errors by iterating
% the results of ps_plot('rsb'). You can add add or subtract an integer number
% of 2pi cycles over a selected region, or correct with respect to the closest
% wrap of another region. You can keep iterating till the 'rsb'
% errors have reduced. If it goes wrong and you do want to reset your 
% unwrapped data back to that of the original you can set the reset_flag to 1;
%
% By Bekaert David - University of Leeds - Sept 2014
% 
% modifications:
% Bekaert David     05/2015     Fix in command line output, fix the code for ps_plot
% Bekaert David     09/2015     Add file restoring.
% Bekaert David     12/2016     Add closest wrap option, change to stamps_save code
% Bekaert David     12/2016     Add argument which is a subset of intergerograms to reduce plotting needs.
%

% deramp IFGs
deramp_flag =0;         % optional deramp rsb. Might be needed when processed with gamma
% what to use as guidence 'rsb' or 'usb' to click on
plot_option = 'usb';    % only 'rsb' or 'usb'
dem_option = 'd';    % only 'd' or ''



% flags
if nargin<1
    reset_flag=0;
end
if isempty(reset_flag)
    reset_flag=0;
end
if nargin<2
    ix_reduced_list=[];
end
if ~strcmpi(plot_option,'usb') && ~strcmpi(plot_option,'rsb')
    error('Only plot_option: rsb and usb supported');
end
if ~strcmpi(dem_option,'d') && ~isempty(dem_option)
    error('Only dem_option: d and [] supported');
end


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


% generating list of ifgs to be ploted
drop_ifg_index = getparm('drop_ifg_index');
if ~isempty(ix_reduced_list)
    % removing ifgs that have been dropped before
    for k=1:length(drop_ifg_index)
        ix_temp = find(drop_ifg_index(k)==ix_reduced_list);
        ix_reduced_list(ix_temp)=[];
    end
    if sum(ix_reduced_list>ps.n_ifg)>0
        fprintf('Your list is larger than number of IFGS, will reset to max number of IFG \n')
        ix_reduced_list(ix_reduced_list>ps.n_ifg)=[];
    end
end
% reset the list in case nothing was left
if isempty(ix_reduced_list)
   ix_reduced_list = 1:ps.n_ifg;
   ix_reduced_list(drop_ifg_index)=[];
end 

% deramping
if deramp_flag==1
    deramp_option = ['-o'];
else
    deramp_option = '';    
end
if ~isempty(dem_option)
    if deramp_flag==1
        dem_option='-do';
    else
        dem_option= '-d';
    end
end
if strcmpi(plot_option,'usb')
    plot_option = [plot_option dem_option];
elseif strcmpi(plot_option,'rsb')
    plot_option = [plot_option deramp_option];
end


% plotting the current rsb data
ps_plot(['rsb' deramp_option],1,0,0,ix_reduced_list);

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
h_fig = ps_plot([plot_option],1,0,0,ix_ifg);
set(h_fig,'name',['Original interferogram']);


% Getting a polygon of the incorrect unwrapped area 
fprintf('Define the incorrect unwrapped region through a polygon by clicking on the figure. \n')
repeat_zoom=1;
while repeat_zoom==1
    action_flag= input('First, zoom to your region of interest. Press [c] to continue. ','s');
    if strcmpi(action_flag,'c') 
        repeat_zoom=0;
    end
end

fprintf('Now, start defining the polygon by outlining the incorrect region. \n Press enter once done\n')
% call the figure in case the user clicked somewhere else
figure(h_fig);
polygon=ginput;
% plotting the polygon on top
hold on
plot([polygon(:,1) ;polygon(1,1)],[polygon(:,2);polygon(1,2)],'r-','linewidth',2)

     
% loop untill the user is happy with it
continue_flag=1;
while continue_flag==1
    repeat=1;
    fprintf('You can shift the whole region by an integer number of cycles or you can put all pixels to a specifc wrap \n')
    while repeat==1
        ix_shift= input('By how many cycles to you want to shift this region? [+-integer or inf for wrap option] ','s');
        ix_shift = str2num(ix_shift);
        if isempty(ix_shift) 
            repeat=1;
        elseif ix_shift==inf
            repeat=0;
        elseif (ix_shift./(round(ix_shift)))~=1
            fprintf(['Needs to be an integer number... \n'])
            repeat =1;
        else
            repeat=0;
        end
    end

    % finding the pixels within the polygon
    ix = inpolygon(ps.lonlat(:,1),ps.lonlat(:,2),polygon(:,1),polygon(:,2));

    % checking the option the user picked - closes wrap (inf) or shift region
    if ix_shift==inf            % closest wrap option
        fprintf('Define region to which you want to define as reference wrap (average will be used!). \n Press enter once done\n')
        % call the figure in case the user clicked somewhere else
        figure(h_fig);
        polygon_ref=ginput;
        % plotting the ploygon on top
        hold on
        plot([polygon_ref(:,1) ;polygon_ref(1,1)],[polygon_ref(:,2);polygon_ref(1,2)],'b-','linewidth',4)

        % finding the pixels within the reference polygon
        ix_ref = inpolygon(ps.lonlat(:,1),ps.lonlat(:,2),polygon_ref(:,1),polygon_ref(:,2));

        % check to which interferograms this should be applied
        repeat2=1;
        while repeat2==1
            action_flag= input('Do you want to apply this to all interferograms [y/n]? ','s');
            if strcmpi(action_flag,'y')
                repeat2=0;
            elseif strcmpi(action_flag,'n')
                repeat2=0;
            end
        end
        
        % store the orginal interferograms
        ix_ifg_or = ix_ifg;
        if strcmpi(action_flag,'y')  
           ix_ifg=[1:size(ph_uw,2)];
        end
        
        % do the estimation for each itnerferogram
        ph_uw= ph_input.ph_uw;
        ref_phase = nanmean(ph_uw(ix_ref,ix_ifg),1);
        for k_ifgs=1:length(ix_ifg)
           % compute the reference 
            radian_shift = round((ph_uw(ix,ix_ifg(k_ifgs))-ref_phase(k_ifgs))./(2*pi))*2*pi;
            ph_uw(ix,ix_ifg(k_ifgs)) = ph_uw(ix,ix_ifg(k_ifgs)) - radian_shift;
        end
        
        % update back to the previous interferogram to be corrected for
        % plotting purposes
        ix_ifg = ix_ifg_or;
        clear ref_phase
        
    else        % option of shifting the interferogram
        % the shift in radians
        radian_shift = ix_shift*2*pi;

        % modifying the interferogram
        ph_uw= ph_input.ph_uw;
        ph_uw(ix,ix_ifg)=ph_uw(ix,ix_ifg)+radian_shift;
    end
    msd = ph_input.msd;

    % saving the data
    stamps_save('phuw_sb2.mat',ph_uw,msd);

    % re-running the 
    sb_invert_uw

    % plot the new residuals
    ps_plot(['rsb' deramp_option],1,0,0,ix_reduced_list);
    h_fig_new = ps_plot([plot_option],1,0,0,ix_ifg);
    set(h_fig_new,'name',['Corrected interferogram']);
    
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
