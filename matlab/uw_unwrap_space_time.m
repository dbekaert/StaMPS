function []=uw_unwrap_space_time(day,unwrap_method,time_win,master_day,bperp,n_trial_wraps)
%UW_UNWRAP_SPACE_TIME smooth and unwrap phase diffs between neighboring data points in time
%
%   Andy Hooper, June 2006
%
%   ============================================================================
%   04/2007 AH: Smoothing changed to time domain (better for non-uniform sampling)
%   11/2009 AH: Extra iteration added on local linear fit for smoothing
%   01/2012 AH: New method 3D_NEW that estimates SCLA for each arc
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

if nargin<5 & strcmpi(unwrap_method,'3D_NEW');
    error('perpendicular baselines must be passed for method 3D_NEW')
end


fprintf('Unwrapping in space-time...\n')

uw=load('uw_grid');
ui=load('uw_interp');

day_pos_ix=find(day>0);
[Y,I]=min(day(day_pos_ix));
close_master_ix=day_pos_ix(I);
if close_master_ix>1
    close_master_ix=[close_master_ix-1;close_master_ix];
end

dph_space=(uw.ph(ui.edges(:,3),:).*conj(uw.ph(ui.edges(:,2),:)));
K=zeros(ui.n_edge,1);

if strcmpi(unwrap_method,'3D_NEW')
    dph_temp=[dph_space(:,[1:close_master_ix-1]),mean(abs(dph_space),2),dph_space(:,[close_master_ix:end])];
    ddph=dph_temp(:,[2:end]).*conj(dph_temp(:,1:end-1)); % sequential dph, to reduce influence of defo
    ddph=ddph./abs(ddph); % normalise
    clear dph_temp
    bperp_master=[bperp(1:close_master_ix-1);0;bperp(close_master_ix:end)];
    bperp_diff=diff(bperp_master);
    bperp_range=max(bperp_diff)-min(bperp_diff);

    trial_mult=[-ceil(8*n_trial_wraps):ceil(8*n_trial_wraps)];
    n_trials=length(trial_mult);
    trial_phase=bperp_diff/bperp_range*pi/4;
    trial_phase_mat=exp(-j*trial_phase*trial_mult);
    coh=zeros(ui.n_edge,1);
    for i=1:ui.n_edge
        cpxphase=ddph(i,:).';
        cpxphase_mat=repmat(cpxphase,1,n_trials);
        phaser=trial_phase_mat.*cpxphase_mat;
        phaser_sum=sum(phaser);
        [coh_max,coh_max_ix]=max(abs(phaser_sum));
        K0=pi/4/bperp_range*trial_mult(coh_max_ix);
        resphase=cpxphase.*exp(-j*(K0*bperp_diff)); % subtract approximate fit
        offset_phase=sum(resphase);
        resphase=angle(resphase*conj(offset_phase)); % subtract offset, take angle (unweighted)
        weighting=abs(cpxphase); 
        mopt=double(weighting.*bperp_diff)\double(weighting.*resphase);
        K(i)=K0+mopt;
        phase_residual=cpxphase.*exp(-j*(K0*bperp_diff)); 
        mean_phase_residual=sum(phase_residual); 
        coh(i)=abs(mean_phase_residual)/sum(abs(phase_residual)); 
    end
    K(coh<0.31)=0;
    dph_space=dph_space.*exp(-j*K*bperp');
end

if strcmpi(unwrap_method,'2D')
    save('uw_space_time','dph_space');
else    
    dph_smooth=zeros(ui.n_edge,uw.n_ifg,'single');
    for i1=1:uw.n_ifg
        time_diff=(day(i1)-day)';
        weight_factor=exp(-(time_diff.^2)/2/time_win^2);
        weight_factor=weight_factor/sum(weight_factor);
            
        dph_mean=sum(dph_space.*repmat(weight_factor,ui.n_edge,1),2);
        dph_mean_adj=angle(dph_space.*repmat(conj(dph_mean),1,uw.n_ifg)); % subtract weighted mean
        G=[ones(uw.n_ifg,1),time_diff'];
        m=lscov(G,double(dph_mean_adj)',weight_factor);
        dph_mean_adj=angle(exp(j*(dph_mean_adj-(G*m)'))); % subtract first estimate
        m2=lscov(G,double(dph_mean_adj)',weight_factor);
        dph_smooth(:,i1)=dph_mean.*exp(j*(m(1,:)'+m2(1,:)')); % add back weighted mean
    end
    dph_noise=angle(dph_space.*conj(dph_smooth));
    dph_smooth_uw=cumsum([angle(dph_smooth(:,1)),angle(dph_smooth(:,2:end).*conj(dph_smooth(:,1:end-1)))],2);

    dph_close_master=mean(dph_smooth_uw(:,close_master_ix),2);
    dph_smooth_uw=dph_smooth_uw-repmat(dph_close_master-angle(exp(j*dph_close_master)),1,uw.n_ifg);
    dph_space_uw=dph_smooth_uw+(dph_noise);
    if strcmpi(unwrap_method,'3D_NEW')
        dph_space_uw=dph_space_uw+K*bperp';
    end


    fprintf('\n   ESTIMATES OF DIFF PHASE NOISE STD DEV\n')
    fprintf('   =====================================\n')
    for i=1:uw.n_ifg
        fprintf('   %s  %4.1f deg\n',datestr(master_day+day(i)),std(dph_noise(:,i))*180/pi)
    end

    save('uw_space_time','dph_space','dph_space_uw','close_master_ix','dph_noise','K');
end
