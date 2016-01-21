function []=uw_unwrap_space_time(day,unwrap_method,time_win,master_day,la_flag,bperp,n_trial_wraps,prefilt_win,scf_flag,temp,n_temp_wraps)
%UW_UNWRAP_SPACE_TIME smooth and unwrap phase diffs between neighboring data points in time
%
%   Andy Hooper, June 2006
%
%   ============================================================================
%   04/2007 AH: Smoothing changed to time domain (better for non-uniform sampling)
%   11/2009 AH: Extra iteration added on local linear fit for smoothing
%   01/2012 AH: New option to estimate look angle error for each arc
%   02/2012 AH: New method 3D_NEW
%   11/2012 AH: Add temperature option
%   ============================================================================

if nargin<1
    help uw_unwrap_time
    error('not enough arguments')
end

if nargin<4
    master_day=0;
end

if nargin<3
    time_win=180;
end

if nargin<2
    unwrap_method='3D';
end

if nargin<5 
    la_flag = 'n';
end

if nargin<6 & strcmpi(la_flag,'y')
    error('perpendicular baselines must be passed for look angle error estimation')
end


fprintf('Unwrapping in space-time...\n')

uw=load('uw_grid');
ui=load('uw_interp');

[nrow,ncol]=size(ui.Z);

day_pos_ix=find(day>0);
[~,I]=min(day(day_pos_ix));
close_master_ix=day_pos_ix(I);
if close_master_ix>1
    close_master_ix=[close_master_ix-1;close_master_ix];
end

dph_space=(uw.ph(ui.edgs(:,3),:).*conj(uw.ph(ui.edgs(:,2),:)));
ifreq_ij=[];
jfreq_ij=[];

if strcmpi(la_flag,'y')
    fprintf('   Estimating look angle error (elapsed time=%ds)\n',round(toc))
    dph_temp=[dph_space(:,[1:close_master_ix-1]),mean(abs(dph_space),2),dph_space(:,[close_master_ix:end])];
    ddph=dph_temp(:,[2:end]).*conj(dph_temp(:,1:end-1)); % sequential dph, to reduce influence of defo
    ddph=ddph./abs(ddph); % normalise
    clear dph_temp
    bperp_master=[bperp(1:close_master_ix-1);0;bperp(close_master_ix:end)];
    bperp_diff=diff(bperp_master);
    bperp_range_orig=max(bperp)-min(bperp);
    bperp_range=max(bperp_diff)-min(bperp_diff);
    n_trial_wraps=n_trial_wraps*(bperp_range_orig/bperp_range);
    ix=bperp_diff~=0;
    bperp_diff=bperp_diff(ix);
    
    trial_mult=[-ceil(8*n_trial_wraps):ceil(8*n_trial_wraps)];
    n_trials=length(trial_mult);
    trial_phase=bperp_diff/bperp_range*pi/4;
    trial_phase_mat=exp(-1i*trial_phase*trial_mult);
    K=zeros(ui.n_edge,1,'single');
    coh=zeros(ui.n_edge,1,'single');
    for i=1:ui.n_edge
        cpxphase=ddph(i,ix).';
        cpxphase_mat=repmat(cpxphase,1,n_trials);
        phaser=trial_phase_mat.*cpxphase_mat;
        phaser_sum=sum(phaser);
        coh_trial=abs(phaser_sum)/sum(abs(cpxphase));
        coh_diff=diff(coh_trial);
        coh_peak_ix=[];
        for i1=2:length(coh_diff)
          if coh_diff(i1)<0 & coh_diff(i1-1)>0
              coh_peak_ix=[coh_peak_ix,i1];
          end
        end 
        coh_peak=coh_trial(coh_peak_ix); % peak values of coherence
        [coh_max,coh_max_peak_ix]=max(coh_peak);
        if length(coh_peak)<2 | min(coh_max-coh_peak([1:coh_max_peak_ix-1,coh_max_peak_ix+1:end]))>0.15 % diff in peaks at least 0.1
            %coh_max_ix=coh_peak_ix(coh_max_peak_ix);
            [~,coh_max_ix]=max(abs(phaser_sum)); % allow for no peaks
            K0=pi/4/bperp_range*trial_mult(coh_max_ix);
            resphase=cpxphase.*exp(-1i*(K0*bperp_diff)); % subtract approximate fit
            offset_phase=sum(resphase);
            resphase=angle(resphase*conj(offset_phase)); % subtract offset, take angle (unweighted)
            weighting=abs(cpxphase); 
            mopt=double(weighting.*bperp_diff)\double(weighting.*resphase);
            K(i)=K0+mopt;
            phase_residual=cpxphase.*exp(-1i*(K(i)*bperp_diff)); 
            mean_phase_residual=sum(phase_residual); 
            coh(i)=abs(mean_phase_residual)/sum(abs(phase_residual)); 
        end
    end
    clear ddph cpxphase_mat trial_phase_mat phaser
    K(coh<0.31)=0;
    dph_space=dph_space.*exp(-1i*K*bperp');   

end

if strcmpi(scf_flag,'y')
    fprintf('   Estimating temperature correlation (elapsed time=%ds)\n',round(toc))
    dph_temp=[dph_space(:,[1:close_master_ix-1]),mean(abs(dph_space),2),dph_space(:,[close_master_ix:end])];
    ddph=dph_temp(:,[2:end]).*conj(dph_temp(:,1:end-1)); % sequential dph, to reduce influence of defo
    ddph=ddph./abs(ddph); % normalise
    clear dph_temp
    temp_master=[temp(1:close_master_ix-1);0;temp(close_master_ix:end)];
    temp_diff=diff(temp_master);
    temp_range_orig=max(temp)-min(temp);
    temp_range=max(temp_diff)-min(temp_diff);
    n_temp_wraps=n_temp_wraps*(temp_range_orig/temp_range);
    ix=temp_diff~=0;
    temp_diff=temp_diff(ix);
    
    trial_mult=[-ceil(8*n_temp_wraps):ceil(8*n_temp_wraps)];
    n_trials=length(trial_mult);
    trial_phase=temp_diff/temp_range*pi/4;
    trial_phase_mat=exp(-1i*trial_phase*trial_mult);
    Kt=zeros(ui.n_edge,1,'single');
    coh=zeros(ui.n_edge,1,'single');
    for i=1:ui.n_edge
        cpxphase=ddph(i,ix).';
        cpxphase_mat=repmat(cpxphase,1,n_trials);
        phaser=trial_phase_mat.*cpxphase_mat;
        phaser_sum=sum(phaser);
        coh_trial=abs(phaser_sum)/sum(abs(cpxphase));
        coh_diff=diff(coh_trial);
        coh_peak_ix=[];
        for i1=2:length(coh_diff)
          if coh_diff(i1)<0 & coh_diff(i1-1)>0
              coh_peak_ix=[coh_peak_ix,i1];
          end
        end 
        coh_peak=coh_trial(coh_peak_ix); % peak values of coherence
        [coh_max,coh_max_peak_ix]=max(coh_peak);
        if length(coh_peak)<2 | min(coh_max-coh_peak([1:coh_max_peak_ix-1,coh_max_peak_ix+1:end]))>0.15 % diff in peaks at least 0.1
            %coh_max_ix=coh_peak_ix(coh_max_peak_ix);
            [~,coh_max_ix]=max(abs(phaser_sum)); % allow for no peaks
            K0=pi/4/temp_range*trial_mult(coh_max_ix);
            resphase=cpxphase.*exp(-1i*(K0*temp_diff)); % subtract approximate fit
            offset_phase=sum(resphase);
            resphase=angle(resphase*conj(offset_phase)); % subtract offset, take angle (unweighted)
            weighting=abs(cpxphase); 
            mopt=double(weighting.*temp_diff)\double(weighting.*resphase);
            Kt(i)=K0+mopt;
            phase_residual=cpxphase.*exp(-1i*(Kt(i)*temp_diff)); 
            mean_phase_residual=sum(phase_residual); 
            coh(i)=abs(mean_phase_residual)/sum(abs(phase_residual)); 
        end
    end
    clear ddph cpxphase_mat trial_phase_mat phaser
    Kt(coh<0.31)=0;
    dph_space=dph_space.*exp(-1i*Kt*temp');   

end

spread=sparse(zeros(ui.n_edge,uw.n_ifg));

if strcmpi(unwrap_method,'2D')
    dph_space_uw=angle(dph_space);
    if strcmpi(la_flag,'y')
        dph_space=dph_space.*exp(1i*K*bperp'); % add back DEM error
        dph_space_uw=dph_space_uw+K*bperp';   % equal to dph_space + integer cycles
    end
    save('uw_space_time','dph_space_uw','spread');    
else
    fprintf('   Smoothing in time (elapsed time=%ds)\n',round(toc))
    dph_smooth=zeros(ui.n_edge,uw.n_ifg,'single');
    dph_space_angle=angle(dph_space);
    for i1=1:uw.n_ifg
        time_diff=(day(i1)-day)';
        weight_factor=exp(-(time_diff.^2)/2/time_win^2);
        weight_factor=weight_factor/sum(weight_factor);
            
        dph_mean=sum(dph_space.*repmat(weight_factor,ui.n_edge,1),2);
        dph_mean_adj=mod(dph_space_angle-repmat(angle(dph_mean),1,uw.n_ifg)+pi,2*pi)-pi;
        G=[ones(uw.n_ifg,1),time_diff'];
        m=lscov(G,double(dph_mean_adj)',weight_factor);
        dph_smooth(:,i1)=dph_mean.*exp(1i*(m(1,:)')); % add back weighted mean
    end
    dph_noise=angle(dph_space.*conj(dph_smooth));
    clear dph_space dph_mean_adj dph_space_angle
    dph_smooth_uw=cumsum([angle(dph_smooth(:,1)),angle(dph_smooth(:,2:end).*conj(dph_smooth(:,1:end-1)))],2);
    clear dph_smooth
    dph_close_master=mean(dph_smooth_uw(:,close_master_ix),2);
    dph_smooth_uw=dph_smooth_uw-repmat(dph_close_master-angle(exp(1i*dph_close_master)),1,uw.n_ifg);
    dph_space_uw=dph_smooth_uw+(dph_noise);
    clear dph_smooth_uw
    
    if strcmpi(la_flag,'y')
        dph_space_uw=dph_space_uw+K*bperp';   % equal to dph_space + integer cycles
    end
    
    if strcmpi(scf_flag,'y')

        fprintf('   Calculating local phase gradients (elapsed time=%ds)\n',round(toc))
        ifreq_ij=nan(uw.n_ps,uw.n_ifg,'single');
        jfreq_ij=nan(uw.n_ps,uw.n_ifg,'single');
        ifgw=zeros(nrow,ncol);

        for i=1:uw.n_ifg
            ifgw(uw.nzix)=uw.ph(:,i);
            [ifreq,jfreq,grad_ij,Hmag]=gradient_filt(ifgw,prefilt_win);
            ix=~isnan(ifreq)&Hmag./(sqrt(abs(ifreq))+1)>3;
            ifreq_ij(:,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),ifreq(ix),uw.ij(:,2),uw.ij(:,1));
            ix=~isnan(jfreq)&Hmag./(sqrt(abs(jfreq))+1)>3;
            jfreq_ij(:,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),jfreq(ix),uw.ij(:,2),uw.ij(:,1));
            nan_ix=isnan(ifreq_ij(:,i));
            ifreq_ij(nan_ix,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),ifreq(ix),uw.ij(nan_ix,2),uw.ij(nan_ix,1),'nearest');
            jfreq_ij(nan_ix,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),jfreq(ix),uw.ij(nan_ix,2),uw.ij(nan_ix,1),'nearest');
        end
        
        fprintf('   Smoothing using local phase gradients (elapsed time=%ds)\n',round(toc))

        dph_smooth_uw2=nan(ui.n_edge,uw.n_ifg,'single');
        spread2=zeros(size(spread),'single');

        for i=1:ui.n_edge
            nodes_ix=ui.edgs(i,[2:3]);
            ifreq_edge=mean(ifreq_ij(nodes_ix,:));
            jfreq_edge=mean(jfreq_ij(nodes_ix,:));
            spread2(i,:)=diff(ifreq_ij(nodes_ix,:))+diff(jfreq_ij(nodes_ix,:));
            dph_smooth_uw2(i,:)=diff(uw.ij(nodes_ix,1))*ifreq_edge+diff(uw.ij(nodes_ix,2))*jfreq_edge;
        end
        fprintf('   Choosing between time and phase gradient smoothing (elapsed time=%ds)\n',round(toc))     
        dph_noise2=angle(exp(1i*(dph_space_uw-dph_smooth_uw2)));
        std_noise=std(dph_noise,0,2);
        std_noise2=std(dph_noise2,0,2);
        dph_noise2(std_noise2>1.3,:)=nan;
        shaky_ix=isnan(std_noise) | std_noise>std_noise2; % spatial smoothing works better index
        fprintf('   %d arcs smoothed in time, %d in space (elapsed time=%ds)\n',ui.n_edge-sum(shaky_ix),sum(shaky_ix),round(toc))        
        dph_noise(shaky_ix,:)=dph_noise2(shaky_ix,:);
        dph_space_uw(shaky_ix,:)=dph_smooth_uw2(shaky_ix,:)+dph_noise2(shaky_ix,:);
        spread(shaky_ix,:)=spread2(shaky_ix,:);
    else
        shaky_ix=[];
    end
    
    
    fprintf('\n   ESTIMATES OF DIFF PHASE NOISE STD DEV\n')
    fprintf('   =====================================\n')
    for i=1:uw.n_ifg
        fprintf('   %s  %4.1f deg\n',datestr(master_day+day(i)),nanstd(dph_noise(:,i))*180/pi)
    end

    dph_noise(std(dph_noise,0,2)>1.3,:)=nan;

    save('uw_space_time','dph_space_uw','close_master_ix','dph_noise','K','spread','ifreq_ij','jfreq_ij','shaky_ix');
end
