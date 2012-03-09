function []=uw_sb_unwrap_space_time(day,ifgday_ix,unwrap_method,time_win,la_flag,bperp,n_trial_wraps,prefilt_win)
%UW_SB_UNWRAP_SPACE_TIME smooth and unwrap phase diffs between neighboring data points in time
%
%   Andy Hooper, June 2007
%
%   ======================================================================
%   07/2008 AH: Allow isolated images, not included in any interferogram
%   01/2010 AH: Allow for non-positive definite inverse var/covar matrix
%   01/2012 AH: New option la_flag that estimates LA error for each arc
%   02/2012 AH: New method 3D_NEW
%   ======================================================================

disp('Unwrapping in time-space...')
tic

uw=load('uw_grid');
ui=load('uw_interp');

n_ifg=uw.n_ifg;
n_ps=uw.n_ps;
n_image=size(day,1);
master_ix=find(day==0);
[nrow,ncol]=size(ui.Z);

day_pos_ix=find(day>0);
[~,I]=min(day(day_pos_ix));
dph_noise=zeros(ui.n_edge,uw.n_ifg,'single');
dph_space=((uw.ph(ui.edges(:,3),:).*conj(uw.ph(ui.edges(:,2),:))));
dph_space=dph_space./abs(dph_space);
K=zeros(ui.n_edge,1);
ifreq_ij=[];
jfreq_ij=[];

if strcmpi(la_flag,'y')
    bperp_range=max(bperp)-min(bperp);

    trial_mult=[-ceil(8*n_trial_wraps):ceil(8*n_trial_wraps)];
    n_trials=length(trial_mult);
    trial_phase=bperp/bperp_range*pi/4;
    trial_phase_mat=exp(-j*trial_phase*trial_mult);
    coh=zeros(ui.n_edge,1);
    for i=1:ui.n_edge
        cpxphase=dph_space(i,:).';
        cpxphase_mat=repmat(cpxphase,1,n_trials);
        phaser=trial_phase_mat.*cpxphase_mat;
        phaser_sum=sum(phaser);
        [~,coh_max_ix]=max(abs(phaser_sum));
        K0=pi/4/bperp_range*trial_mult(coh_max_ix);
        resphase=cpxphase.*exp(-1i*(K0*bperp)); % subtract approximate fit
        offset_phase=sum(resphase);
        resphase=angle(resphase*conj(offset_phase)); % subtract offset, take angle (unweighted)
        weighting=abs(cpxphase); 
        mopt=double(weighting.*bperp)\double(weighting.*resphase);
        K(i)=K0+mopt;
        phase_residual=cpxphase.*exp(-1i*(K0*bperp)); 
        mean_phase_residual=sum(phase_residual); 
        coh(i)=abs(mean_phase_residual)/sum(abs(phase_residual)); 
    end
    dph_space=dph_space.*exp(-1i*K*bperp');
end

spread=zeros(ui.n_edge,n_ifg);

if strcmpi(unwrap_method,'2D')
    dph_space_uw=angle(dph_space);
    if strcmpi(la_flag,'y')
        dph_space=dph_space.*exp(1i*K*bperp'); % add back DEM error
        dph_space_uw=dph_space_uw+K*bperp';   % equal to dph_space + integer cycles
    end
    save('uw_space_time','dph_space','dph_space_uw','spread');    
elseif strcmpi(unwrap_method,'3D_NO_DEF')
    dph_noise=angle(dph_space);
    dph_space_uw=angle(dph_space);        
    if strcmpi(la_flag,'y')
        dph_space=dph_space.*exp(1i*K*bperp'); % add back DEM error
        dph_space_uw=dph_space_uw+K*bperp';   % equal to dph_space + integer cycles
    end
    save('uw_space_time','dph_space','dph_space_uw','dph_noise','spread');
else

    G=zeros(n_ifg,n_image);
    for i=1:n_ifg
        G(i,ifgday_ix(i,1))=-1;
        G(i,ifgday_ix(i,2))=1;
    end

    nzc_ix=sum(abs(G))~=0; % non-zero column index
    day=day(nzc_ix);
    G=G(:,nzc_ix);
    n=size(G,2);
    
    %x=(0:n-1)'; % use sequence for smoothing
    x=(day-day(1))*(n-1)/(day(end)-day(1)); % use dates for smoothing
    
    dph_space_series=[zeros(1,ui.n_edge);double(G(:,2:end))\double(angle(dph_space))'];
    dph_smooth_series=zeros(size(G,2),ui.n_edge,'single');

    for i1=1:n
        time_diff_sq=(day(i1)-day).^2;
        weight_factor=exp(-time_diff_sq/2/time_win^2);
        weight_factor=weight_factor/sum(weight_factor);
        dph_smooth_series(i1,:)=sum(dph_space_series.*repmat(weight_factor,1,ui.n_edge));
    end

    dph_smooth_ifg=(G*dph_smooth_series)';
    dph_noise=angle(dph_space.*exp(-1i*dph_smooth_ifg));

    if strcmpi(unwrap_method,'3D_SMALL_DEF')|...
       strcmpi(unwrap_method,'3D_QUICK')|...
       strcmpi(unwrap_method,'3D_NEW')
      not_small_ix=find(std(dph_noise,0,2)>1.3)';
      fprintf('Ignoring %d edges. Elapsed time = %d s\n',length(not_small_ix),round(toc))
      dph_noise(not_small_ix,:)=nan;
    else
      ph_noise=angle(uw.ph.*conj(uw.ph_lowpass));
      dph_noise_sf=((ph_noise(ui.edges(:,3),:)-(ph_noise(ui.edges(:,2),:))));      
      m_minmax=repmat([-pi,pi],5,1).*repmat([0.5;0.25;1;0.25;1],1,2);
      anneal_opts=[1;15;0;0;0;0;0];
      covm=cov((dph_noise_sf)); % estimate of covariance
      [W,P]=chol(inv(covm)); % weighting matrix 
      if P~=0
          W=diag(1./sqrt(diag(covm)));
      end
      not_small_ix=find(std(dph_noise,0,2)>1)';
      fprintf('Performing more complex smoothing on %d edges. Elapsed time = %d s\n',length(not_small_ix),round(toc))

      n_proc=0;
      for i=not_small_ix
        dph=angle(dph_space(i,:))';
        dph_smooth_series(:,i)=uw_sb_smooth_unwrap(m_minmax,anneal_opts,G,W,dph,x);

        n_proc=n_proc+1;
        if round(n_proc/1000)==n_proc/1000
            save('uw_unwrap_time','G','dph_space','dph_smooth_series');
            fprintf('%d edges of %d reprocessed. Elapsed time = %d s\n',n_proc,length(not_small_ix),round(toc))
        end
      end
      dph_smooth_ifg=(G*dph_smooth_series)';
      dph_noise=angle(dph_space.*exp(-1i*dph_smooth_ifg));
    end
    dph_space_uw=dph_smooth_ifg+dph_noise;
    if strcmpi(la_flag,'y')
        dph_space=dph_space.*exp(1i*K*bperp'); % add back DEM error
        dph_space_uw=dph_space_uw+K*bperp';   % equal to dph_space + integer cycles
    end

    if strcmpi(unwrap_method,'3D_NEW')

        ifreq_ij=nan(n_ps,n_ifg);
        jfreq_ij=nan(n_ps,n_ifg);
        ifgw=zeros(nrow,ncol);
        dph_smooth_uw2=nan(ui.n_edge,n_ifg);
        spread2=spread;
        for i=1:n_ifg
            ifgw(uw.nzix)=uw.ph(:,i);
            [ifreq,jfreq,grad_ij,Hmag]=gradient_filt(ifgw,prefilt_win);
            ix=~isnan(ifreq)&Hmag>3;
            ifreq_ij(:,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),ifreq(ix),uw.ij(:,2),uw.ij(:,1));
            jfreq_ij(:,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),jfreq(ix),uw.ij(:,2),uw.ij(:,1));
            nan_ix=isnan(ifreq_ij(:,i));
            ifreq_ij(nan_ix,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),ifreq(ix),uw.ij(nan_ix,2),uw.ij(nan_ix,1),'nearest');
            jfreq_ij(nan_ix,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),jfreq(ix),uw.ij(nan_ix,2),uw.ij(nan_ix,1),'nearest');
        end

        for i=1:ui.n_edge
            nodes_ix=ui.edges(i,[2:3]);
            ifreq_edge=mean(ifreq_ij(nodes_ix,:));
            jfreq_edge=mean(jfreq_ij(nodes_ix,:));
            spread2(i,:)=diff(ifreq_ij(nodes_ix,:))+diff(jfreq_ij(nodes_ix,:));
            dph_smooth_uw2(i,:)=diff(uw.ij(nodes_ix,1))*ifreq_edge+diff(uw.ij(nodes_ix,2))*jfreq_edge;
        end
        dph_noise2=angle((dph_space).*exp(-j*dph_smooth_uw2));
        std_noise=std(dph_noise,0,2);
        std_noise2=std(dph_noise2,0,2);
        dph_noise2(std_noise2>1.3,:)=nan;
        shaky_ix=isnan(std_noise) | std_noise>std_noise2; % spatial smoothing works better index
        shaky_nodes=ui.edges(shaky_ix,[2:3]);
        shaky_nodes=sort(shaky_nodes(:));
        not_uniq_ix=find(diff(shaky_nodes)==0);
        shaky_nodes=shaky_nodes(not_uniq_ix);
        shaky_nodes=unique(shaky_nodes);
        for i=1:length(shaky_nodes)
            shaky_edges=(ui.edges(:,2)==shaky_nodes(i)|ui.edges(:,3)==shaky_nodes(i));
            shaky_ix(shaky_edges)=true; % for nodes with >1 shaky edges, set all edges to shaky
        end

        dph_noise(shaky_ix,:)=dph_noise2(shaky_ix,:);
        dph_space_uw(shaky_ix,:)=dph_smooth_uw2(shaky_ix,:)+dph_noise2(shaky_ix,:);
        spread(shaky_ix,:)=spread2(shaky_ix,:);
    end
    
    save('uw_space_time','dph_space','dph_space_uw','dph_noise','G','dph_smooth_series','spread','ifreq_ij','jfreq_ij');
end
