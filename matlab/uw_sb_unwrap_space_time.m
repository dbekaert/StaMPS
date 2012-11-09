function []=uw_sb_unwrap_space_time(day,ifgday_ix,unwrap_method,time_win,la_flag,bperp,n_trial_wraps,prefilt_win,scf_flag,temp,n_temp_wraps)
%UW_SB_UNWRAP_SPACE_TIME smooth and unwrap phase diffs between neighboring data points in time
%
%   Andy Hooper, June 2007
%
%   ======================================================================
%   07/2008 AH: Allow isolated images, not included in any interferogram
%   01/2010 AH: Allow for non-positive definite inverse var/covar matrix
%   01/2012 AH: New option la_flag that estimates LA error for each arc
%   02/2012 AH: New method 3D_NEW
%   02/2012 AH: New method 3D_FULL
%   03/2012 AH: version updated for subsets of full network
%   10/2012 AH: Bug corrected in 3D_FULL
%   11/2012 AH: Add temperature option and scf_flag
%   ======================================================================
% 
%
disp('Unwrapping in time-space...')


uw=load('uw_grid');
ui=load('uw_interp');

n_ifg=uw.n_ifg;
n_ps=uw.n_ps;
nzix=uw.nzix;
ij=uw.ij;

n_image=size(day,1);
master_ix=find(day==0);
[nrow,ncol]=size(ui.Z);

day_pos_ix=find(day>0);
[~,I]=min(day(day_pos_ix));
dph_space=((uw.ph(ui.edges(:,3),:).*conj(uw.ph(ui.edges(:,2),:))));
clear uw

dph_space=dph_space./abs(dph_space);
ifreq_ij=[];
jfreq_ij=[];

G=zeros(n_ifg,n_image);
for i=1:n_ifg
    G(i,ifgday_ix(i,1))=-1;
    G(i,ifgday_ix(i,2))=1;
end

nzc_ix=sum(abs(G))~=0; % non-zero column index
day=day(nzc_ix);
n_image=size(day,1);
G=G(:,nzc_ix);
zc_ix=find(nzc_ix==0);
for i=1:length(zc_ix)
    ifgday_ix(ifgday_ix>zc_ix(i))=ifgday_ix(ifgday_ix>zc_ix(i))-1;
end
n=size(G,2);

if ~isempty(temp)
    temp_flag='y';
else
    temp_flag='n';
end
    
if strcmpi(la_flag,'y') | strcmpi(temp_flag,'y')
    
    if isempty(bperp)
        bperp=[1:n_ifg]';
    end
    if isempty(temp)
        temp=[1:n_ifg]';
    end
    
    bperp_range=max(bperp)-min(bperp);
    temp_range=max(temp)-min(temp);
    ix=find(abs(diff(ifgday_ix,[],2))==1);   

    if length(ix)>=length(day)-1 % if almost full cascade is present
        fprintf('   using sequential cascade of interferograms\n')
        temp_sub=temp(ix);
        temp_range_sub=max(temp_sub)-min(temp_sub);
        dph_sub=dph_space(:,ix); % only ifgs using ith image
        bperp_sub=bperp(ix);
        bperp_range_sub=max(bperp_sub)-min(bperp_sub);
        n_trial_wraps=n_trial_wraps*(bperp_range_sub/bperp_range);
        n_temp_wraps=n_temp_wraps*(temp_range_sub/temp_range);
    else
        ifgs_per_image=sum(abs(G));
        [max_ifgs_per_image,max_ix]=max(ifgs_per_image);
        if max_ifgs_per_image>=length(day)-2; % can create cascade from (almost) complete single master series
            fprintf('   using sequential cascade of interferograms\n')
            ix=G(:,max_ix)~=0;
            gsub=G(ix,max_ix);
            sign_ix=-sign(single(gsub'));
            dph_sub=dph_space(:,ix); % only ifgs using chosen master image
            bperp_sub=[bperp(ix)]; 
            bperp_sub(sign_ix==-1)=-bperp_sub(sign_ix==-1);
            bperp_sub=[bperp_sub;0]; % add master
            temp_sub=[temp(ix)]; 
            temp_sub(sign_ix==-1)=-temp_sub(sign_ix==-1);
            temp_sub=[temp_sub;0]; % add master
            sign_ix=repmat(sign_ix,ui.n_edge,1);
            dph_sub(sign_ix==-1)=conj(dph_sub(sign_ix==-1)); % flip sign if necessary to make ith image master 
            dph_sub=[dph_sub,mean(abs(dph_sub),2)]; % add zero phase master
            slave_ix=sum(ifgday_ix(ix,:),2)-max_ix;
            day_sub=day([slave_ix;max_ix]); % extract days for subset
            [day_sub,sort_ix]=sort(day_sub); % sort ascending day
            dph_sub=dph_sub(:,sort_ix); % sort ascending day
            bperp_sub=bperp_sub(sort_ix);
            bperp_sub=diff(bperp_sub);
            bperp_range_sub=max(bperp_sub)-min(bperp_sub);
            n_trial_wraps=n_trial_wraps*(bperp_range_sub/bperp_range);           
            temp_sub=temp_sub(sort_ix);
            temp_sub=diff(temp_sub);
            temp_range_sub=max(temp_sub)-min(temp_sub);
            n_temp_wraps=n_temp_wraps*(temp_range_sub/temp_range);           
            n_sub=length(day_sub);
            dph_sub=dph_sub(:,[2:end]).*conj(dph_sub(:,1:end-1)); % sequential dph, to reduce influence of defo
            dph_sub=dph_sub./abs(dph_sub); % normalise
        else % just use all available
            dph_sub=dph_space;
            bperp_sub=bperp;
            bperp_range_sub=bperp_range;
            temp_sub=temp;
            temp_range_sub=temp_range;
        end
    end
end
    
if strcmpi(la_flag,'y') 
    fprintf('   Estimating look angle error (elapsed time=%ds)\n',round(toc))
    
    trial_mult=[-ceil(8*n_trial_wraps):ceil(8*n_trial_wraps)];
    n_trials=length(trial_mult);
    trial_phase=bperp_sub/bperp_range_sub*pi/4;
    trial_phase_mat=exp(-j*trial_phase*trial_mult);
    K=zeros(ui.n_edge,1,'single');
    coh=zeros(ui.n_edge,1,'single');
    for i=1:ui.n_edge
        cpxphase=dph_sub(i,:).';
        cpxphase_mat=repmat(cpxphase,1,n_trials);
        phaser=trial_phase_mat.*cpxphase_mat;
        phaser_sum=sum(phaser);
        coh_trial=abs(phaser_sum)/sum(abs(cpxphase));
        [coh_max,coh_max_ix]=max(coh_trial); % allow for no peaks
        peak_ix=coh_max_ix-8:coh_max_ix+8;
        peak_ix=peak_ix(peak_ix>0&peak_ix<=n_trials);
        coh_trial(peak_ix)=0;
        if coh_max-max(coh_trial)>0.1 % diff in peaks at least 0.1
%         coh_diff=diff(coh_trial);
%         coh_peak_ix=[];
%         for i1=2:length(coh_diff)
%           if coh_diff(i1)<0 & coh_diff(i1-1)>0
%               coh_peak_ix=[coh_peak_ix,i1];
%           end
%         end 
%         coh_peak=coh_trial(coh_peak_ix); % peak values of coherence
%         [coh_max,coh_max_peak_ix]=max(coh_peak);
%         if length(coh_peak)<2 | min(coh_max-coh_peak([1:coh_max_peak_ix-1,coh_max_peak_ix+1:end]))>0.1 % diff in peaks at least 0.1
            %coh_max_ix=coh_peak_ix(coh_max_peak_ix);
%            [~,coh_max_ix]=max(abs(phaser_sum)); % allow for no peaks
            K0=pi/4/bperp_range_sub*trial_mult(coh_max_ix);
            resphase=cpxphase.*exp(-1i*(K0*bperp_sub)); % subtract approximate fit
            offset_phase=sum(resphase);
            resphase=angle(resphase*conj(offset_phase)); % subtract offset, take angle (unweighted)
            weighting=abs(cpxphase); 
            mopt=double(weighting.*bperp_sub)\double(weighting.*resphase);
            K(i)=K0+mopt;
            phase_residual=cpxphase.*exp(-1i*(K(i)*bperp_sub)); 
            mean_phase_residual=sum(phase_residual); 
            coh(i)=abs(mean_phase_residual)/sum(abs(phase_residual)); 
        end
    end
    
    clear cpxphase_mat trial_phase_mat phaser 
    K(coh<0.31)=0; % not to be trusted;
    dph_space=dph_space.*exp(-1i*K*bperp');
    dph_sub=dph_sub.*exp(-1i*K*bperp_sub');
        
end


if strcmpi(temp_flag,'y')
    fprintf('   Estimating temperature correlation (elapsed time=%ds)\n',round(toc))

    trial_mult=[-ceil(8*n_temp_wraps):ceil(8*n_temp_wraps)];
    n_trials=length(trial_mult);
    trial_phase=temp_sub/temp_range_sub*pi/4;
    trial_phase_mat=exp(-j*trial_phase*trial_mult);
    Kt=zeros(ui.n_edge,1,'single');
    coh=zeros(ui.n_edge,1,'single');
    for i=1:ui.n_edge
        cpxphase=dph_sub(i,:).';
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
        if length(coh_peak)<2 | min(coh_max-coh_peak([1:coh_max_peak_ix-1,coh_max_peak_ix+1:end]))>0.1 % diff in peaks at least 0.1
            [~,coh_max_ix]=max(abs(phaser_sum)); % allow for no peaks
            K0=pi/4/temp_range_sub*trial_mult(coh_max_ix);
            resphase=cpxphase.*exp(-1i*(K0*temp_sub)); % subtract approximate fit
            offset_phase=sum(resphase);
            resphase=angle(resphase*conj(offset_phase)); % subtract offset, take angle (unweighted)
            weighting=abs(cpxphase); 
            mopt=double(weighting.*temp_sub)\double(weighting.*resphase);
            Kt(i)=K0+mopt;
            phase_residual=cpxphase.*exp(-1i*(Kt(i)*temp_sub)); 
            mean_phase_residual=sum(phase_residual); 
            coh(i)=abs(mean_phase_residual)/sum(abs(phase_residual)); 
        end
    end
        
    clear cpxphase_mat trial_phase_mat phaser dph_sub
    Kt(coh<0.31)=0; % not to be trusted;
    dph_space=dph_space.*exp(-1i*Kt*temp');
        
end

spread=sparse(zeros(ui.n_edge,n_ifg));

if strcmpi(unwrap_method,'2D')
    dph_space_uw=angle(dph_space);
    if strcmpi(la_flag,'y')
        dph_space_uw=dph_space_uw+K*bperp';   % equal to dph_space + integer cycles
    end
    if ~isempty(temp)
        dph_space_uw=dph_space_uw+Kt*temp';   % equal to dph_space + integer cycles
    end
    save('uw_space_time','dph_space_uw','spread');    
elseif strcmpi(unwrap_method,'3D_NO_DEF')
    dph_noise=angle(dph_space);
    dph_space_uw=angle(dph_space);        
    if strcmpi(la_flag,'y')
        dph_space_uw=dph_space_uw+K*bperp';   % equal to dph_space + integer cycles
    end
    if ~isempty(temp)
        dph_space_uw=dph_space_uw+Kt*temp';   % equal to dph_space + integer cycles
    end
    save('uw_space_time','dph_space_uw','dph_noise','spread');
else
    fprintf('   Smoothing in time (elapsed time=%ds)\n',round(toc))
    
    if strcmpi(unwrap_method,'3D_FULL') 
        
       dph_smooth_ifg=nan(size(dph_space),'single');
       
       for i=1:n_image
         ix=G(:,i)~=0;
         if sum(ix)>=n_image-2 % allow only 1 missing ifg for each image as master
           gsub=G(ix,i);
           dph_sub=dph_space(:,ix); % only ifgs using ith image
           sign_ix=repmat(-sign(single(gsub')),ui.n_edge,1);
           dph_sub(sign_ix==-1)=conj(dph_sub(sign_ix==-1)); % flip sign if necessary to make ith image master 
           slave_ix=sum(ifgday_ix(ix,:),2)-i;
           day_sub=day(slave_ix); % extract days for subset
           [day_sub,sort_ix]=sort(day_sub); % sort ascending day
           dph_sub=dph_sub(:,sort_ix); % sort ascending day
           dph_sub_angle=angle(dph_sub);
           n_sub=length(day_sub);
           dph_smooth=zeros(ui.n_edge,n_sub,'single');
           for i1=1:n_sub
               time_diff=(day_sub(i1)-day_sub)';
               weight_factor=exp(-(time_diff.^2)/2/time_win^2);
               weight_factor=weight_factor/sum(weight_factor);

               dph_mean=sum(dph_sub.*repmat(weight_factor,ui.n_edge,1),2);
               %dph_mean_adj=angle(dph_sub.*repmat(conj(dph_mean),1,n_sub)); % subtract weighted mean
               dph_mean_adj=mod(dph_sub_angle-repmat(angle(dph_mean),1,n_sub)+pi,2*pi)-pi;
               GG=[ones(n_sub,1),time_diff'];
               if size(GG,1)>1
                   m=lscov(GG,double(dph_mean_adj)',weight_factor);
               else
                   m=zeros(size(GG,1),ui.n_edge);
               end
               %dph_mean_adj=mod(dph_mean_adj-(GG*m)'+pi,2*pi)-pi; % subtract first estimate
               %m2=lscov(GG,double(dph_mean_adj)',weight_factor);
               %dph_smooth(:,i1)=dph_mean.*exp(1i*(m(1,:)'+m2(1,:)')); % add back weighted mean
               dph_smooth(:,i1)=dph_mean.*exp(1i*(m(1,:)')); % add back weighted mean
            end
            dph_smooth_sub=cumsum([angle(dph_smooth(:,1)),angle(dph_smooth(:,2:end).*conj(dph_smooth(:,1:end-1)))],2);
            close_master_ix=find(slave_ix-i>0);
            if isempty(close_master_ix)
                close_master_ix=n_sub;
            else
                close_master_ix=close_master_ix(1);
                if close_master_ix>1
                    close_master_ix=[close_master_ix-1;close_master_ix];
                end
            end
            dph_close_master=mean(dph_smooth_sub(:,close_master_ix),2);
            dph_smooth_sub=dph_smooth_sub-repmat(dph_close_master-angle(exp(j*dph_close_master)),1,n_sub);
            dph_smooth_sub=dph_smooth_sub.*sign_ix;
            already_sub_ix=find(~isnan(dph_smooth_ifg(1,ix))); % already unwrapped index
            ix=find(ix);
            already_ix=ix(already_sub_ix);
            std_noise1=std(angle(dph_space(:,already_ix).*exp(-1i*dph_smooth_ifg(:,already_ix))));
            std_noise2=std(angle(dph_space(:,already_ix).*exp(-1i*dph_smooth_sub(:,already_sub_ix))));
            keep_ix=true(n_sub,1);
            keep_ix(already_sub_ix(std_noise1<std_noise2))=false; % keep least noisy
            %keep_ix(already_sub_ix(std_noise1>std_noise2))=false; % keep most noisy
            dph_smooth_ifg(:,ix(keep_ix))=dph_smooth_sub(:,keep_ix);
          end
        end
       
        dph_noise=angle(dph_space.*exp(-1i*dph_smooth_ifg));
        dph_noise(std(dph_noise,0,2)>1.2,:)=nan;
     
    else
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
           strcmpi(unwrap_method,'3D_QUICK')
          not_small_ix=find(std(dph_noise,0,2)>1.3)';
          fprintf('   %d edges with high std dev in time (elapsed time=%ds)\n',length(not_small_ix),round(toc))
          dph_noise(not_small_ix,:)=nan;
        else % 3D
          uw=load('uw_grid');
          ph_noise=angle(uw.ph.*conj(uw.ph_lowpass));
          clear uw
          dph_noise_sf=((ph_noise(ui.edges(:,3),:)-(ph_noise(ui.edges(:,2),:))));      
          m_minmax=repmat([-pi,pi],5,1).*repmat([0.5;0.25;1;0.25;1],1,2);
          anneal_opts=[1;15;0;0;0;0;0];
          covm=cov((dph_noise_sf)); % estimate of covariance
          [W,P]=chol(inv(covm)); % weighting matrix 
          if P~=0
              W=diag(1./sqrt(diag(covm)));
          end
          not_small_ix=find(std(dph_noise,0,2)>1)';
          fprintf('   Performing complex smoothing on %d edges (elapsed time=%ds)\n',length(not_small_ix),round(toc))

          n_proc=0;
          for i=not_small_ix
            dph=angle(dph_space(i,:))';
            dph_smooth_series(:,i)=uw_sb_smooth_unwrap(m_minmax,anneal_opts,G,W,dph,x);

            n_proc=n_proc+1;
            if round(n_proc/1000)==n_proc/1000
                save('uw_unwrap_time','G','dph_space','dph_smooth_series');
                fprintf('%d edges of %d reprocessed (elapsed time=%ds)\n',n_proc,length(not_small_ix),round(toc))
            end
          end
          dph_smooth_ifg=(G*dph_smooth_series)';
          dph_noise=angle(dph_space.*exp(-1i*dph_smooth_ifg));
        end
    end
    clear dph_space
    dph_space_uw=dph_smooth_ifg+dph_noise;
    clear dph_smooth_ifg
    
    if strcmpi(la_flag,'y')
        dph_space_uw=dph_space_uw+K*bperp';   % equal to dph_space + integer cycles
    end
    if strcmpi(temp_flag,'y')
        dph_space_uw=dph_space_uw+Kt*temp';   % equal to dph_space + integer cycles
    end
    
    if strcmpi(scf_flag,'y')
       
        fprintf('   Calculating local phase gradients (elapsed time=%ds)\n',round(toc))
        ifreq_ij=nan(n_ps,n_ifg,'single');
        jfreq_ij=nan(n_ps,n_ifg,'single');
        ifgw=zeros(nrow,ncol);
        uw=load('uw_grid');
        for i=1:n_ifg
            ifgw(nzix)=uw.ph(:,i);
            [ifreq,jfreq,grad_ij,Hmag]=gradient_filt(ifgw,prefilt_win);
            ix=~isnan(ifreq)&Hmag./((abs(ifreq))+1)>3;
            if sum(ix(:))>2
                ifreq_ij(:,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),ifreq(ix),ij(:,2),ij(:,1),'linear');
            end
            ix=~isnan(jfreq)&Hmag./((abs(jfreq))+1)>3;
            if sum(ix(:))>2
                jfreq_ij(:,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),jfreq(ix),ij(:,2),ij(:,1),'linear');
            end
            %nan_ix=isnan(ifreq_ij(:,i));
            %ifreq_ij(nan_ix,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),ifreq(ix),ij(nan_ix,2),ij(nan_ix,1),'nearest');
            %jfreq_ij(nan_ix,i)=griddata(grad_ij(ix,2),grad_ij(ix,1),jfreq(ix),ij(nan_ix,2),ij(nan_ix,1),'nearest');
        end
        clear uw
        
        spread2=zeros(size(spread),'single');
        dph_smooth_uw2=nan(ui.n_edge,n_ifg,'single');
       
        fprintf('   Smoothing using local phase gradients (elapsed time=%ds)\n',round(toc))
        for i=1:ui.n_edge
            nodes_ix=ui.edges(i,[2:3]);
            ifreq_edge=mean(ifreq_ij(nodes_ix,:));
            jfreq_edge=mean(jfreq_ij(nodes_ix,:));
            diff_i=diff(ij(nodes_ix,1));
            diff_j=diff(ij(nodes_ix,2));
            dph_smooth_uw2(i,:)=diff_i*ifreq_edge+diff_j*jfreq_edge;
%           spread2(i,:)=diff_i*diff(ifreq_ij(nodes_ix,:))+diff_j*diff(jfreq_ij(nodes_ix,:));
            spread2(i,:)=diff(ifreq_ij(nodes_ix,:))+diff(jfreq_ij(nodes_ix,:));
        end
        fprintf('   Choosing between time and phase gradient smoothing (elapsed time=%ds)\n',round(toc))        
        std_noise=std(dph_noise,0,2);
        dph_noise2=angle(exp(j*(dph_space_uw-dph_smooth_uw2)));
        std_noise2=std(dph_noise2,0,2);
        dph_noise2(std_noise2>1.3,:)=nan;
        shaky_ix=isnan(std_noise) | std_noise>std_noise2; % spatial smoothing works better index
        
        %%%FIX FIX
        %shaky_ix=isnan(std_noise) | std_noise>0; % spatial smoothing always works better index
        %shaky_ix= std_noise<0; % spatial smoothing never works better index
        %%%%%%%%%%
        
        %shaky_nodes=ui.edges(shaky_ix,[2:3]);
        %shaky_nodes=sort(shaky_nodes(:));
        %not_uniq_ix= diff(shaky_nodes)==0;
        %shaky_nodes=shaky_nodes(not_uniq_ix);
        %shaky_nodes=unique(shaky_nodes);
        %for i=1:leph_uw=uw_3d(ph(:,ix),xy,day,ifgday_ix(ix,:),bperp(ix),options);ngth(shaky_nodes)
        %    shaky_edges=(ui.edges(:,2)==shaky_nodes(i)|ui.edges(:,3)==shaky_nodes(i));
        %    shaky_ix(shaky_edges)=true; % for nodes with >1 shaky edges, set all edges to shaky
        %end
        fprintf('   %d arcs smoothed in time, %d in space (elapsed time=%ds)\n',ui.n_edge-sum(shaky_ix),sum(shaky_ix),round(toc))        

        dph_noise(shaky_ix,:)=dph_noise2(shaky_ix,:);
        dph_space_uw(shaky_ix,:)=dph_smooth_uw2(shaky_ix,:)+dph_noise2(shaky_ix,:);
        spread(shaky_ix,:)=spread2(shaky_ix,:);
    else
        shaky_ix=[];
    end
    
    save('uw_space_time','dph_space_uw','dph_noise','G','spread','ifreq_ij','jfreq_ij','shaky_ix');
end
