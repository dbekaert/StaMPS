function []=uw_unwrap_space_time(day,unwrap_method,time_win,master_day)
%UW_UNWRAP_SPACE_TIME smooth and unwrap phase diffs between neighboring data points in time
%
%   Andy Hooper, June 2006
%
%   ============================================================================
%   04/2007 AH: Smoothing changed to time domain (better for non-uniform sampling)
%   11/2009 AH: Extra iteration added on local linear fit for smoothing
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
        WG=G.*[weight_factor',weight_factor'];
        m=double(WG)\(repmat(double(weight_factor),ui.n_edge,1).*double(dph_mean_adj))'; % weighted least-sq to find best-fit line
        dph_mean_adj=angle(exp(j*(dph_mean_adj-(G*m)'))); % subtract first estimate
        m2=double(WG)\(repmat(double(weight_factor),ui.n_edge,1).*double(dph_mean_adj))'; % weighted least-sq to find best-fit line
        dph_smooth(:,i1)=dph_mean.*exp(j*(m(1,:)'+m2(1,:)')); % add back weighted mean
        
    end
    dph_noise=angle(dph_space.*conj(dph_smooth));
    dph_smooth_uw=cumsum([angle(dph_smooth(:,1)),angle(dph_smooth(:,2:end).*conj(dph_smooth(:,1:end-1)))],2);

    dph_close_master=mean(dph_smooth_uw(:,close_master_ix),2);
    dph_smooth_uw=dph_smooth_uw-repmat(dph_close_master-angle(exp(j*dph_close_master)),1,uw.n_ifg);
    dph_space_uw=dph_smooth_uw+(dph_noise);


    fprintf('\n   ESTIMATES OF DIFF PHASE NOISE STD DEV\n')
    fprintf('   =====================================\n')
    for i=1:uw.n_ifg
        fprintf('   %s  %4.1f deg\n',datestr(master_day+day(i)),std(dph_noise(:,i))*180/pi)
    end

    save('uw_space_time','dph_space','dph_space_uw','close_master_ix','dph_noise');
end
