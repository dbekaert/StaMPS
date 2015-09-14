function []=sb_identify_good_pixels()
%SB_IDENTIFY_GOOD_PIXELS
%   Do the closure tests to identify 'good pixels', i.e. ones that are
%   likely unwrapped properly.
%   Ekbal Hussain @ Uni Leeds
%   19/11/2013
%
%   ======================================================================
%   06/2014 AH: Modify for StaMPS and fix
%   04/2015 EH: Speed optimisations
%   09/2015 EH/DB: Account for Nan values
%   ======================================================================
logit;
fprintf('Identifying good pixels...\n')

load psver
psname=['ps',num2str(psver)];
phuwname=['phuw_sb',num2str(psver)];
loopname=['phuw_loops',num2str(psver)];
goodname=['phuw_good',num2str(psver)];

sb_make_closure_loops;

ps=load(psname);
load(loopname);
load(phuwname,'ph_uw')

good_pixels = false(ps.n_ps,ps.n_ifg);


for m=1:size(intfg_loops,1)
    
    %Sum positive interferograms in each loop
    positive = find(intfg_loops(m,:)==1);   
    positive_ints = zeros(size(ph_uw,1),1);
    for p=1:size(positive,2)
       positive_ints = positive_ints + ph_uw(:,positive(1,p));
    end
    
    %Sum negative interferograms in each loop    
    negative = find(intfg_loops(m,:)==-1);
    negative_ints = zeros(size(ph_uw,1),1);
    for p=1:size(negative,2)
       negative_ints = negative_ints + ph_uw(:,negative(1,p));
    end
    
    ints_used = [positive, negative];
    
    %Remove the average residual (closest 2*pi value) to ensure each int is
    %in the same reference
    average_resid = nanmean(positive_ints-negative_ints);

    average_resid = 2*pi*(round(average_resid/2/pi));
   
    phase_resid = positive_ints-negative_ints-average_resid;    
    %hist(positive_ints-negative_ints,100)
    %a = phase_resid(find(phase_resid(:,1)<=1 & phase_resid(:,1)>=-1),:);
    %size(a,1)
    
    %Index each interferogram with the pixels that are good=1 and bad=0
    %A good pixel is good if: -1<closure residual<1
    good_pixels(phase_resid(:,1)<=1 & phase_resid(:,1)>=-1,ints_used)=true;      
end

save(goodname,'good_pixels')
logit(1);
