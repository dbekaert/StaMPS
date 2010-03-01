function []=uw_sb_unwrap_time(day,ifgday_ix,unwrap_method,time_win)
%UW_SB_UNWRAP_SPACE_TIME smooth and unwrap phase diffs between neighboring data points in time
%
%   Andy Hooper, June 2007
%
%   ======================================================================
%   07/2008 AH: Allow isolated images, not included in any interferogram
%   01/2010 AH: Allow for non-positive definite inverse var/covar matrix
%   ======================================================================

tic

disp('Unwrapping in time-space...')

uw=load('uw_grid');
ui=load('uw_interp');

n_ifg=uw.n_ifg;
n_image=size(day,1);
master_ix=find(day==0);

day_pos_ix=find(day>0);
[Y,I]=min(day(day_pos_ix));
%close_master_ix=day_pos_ix(I);
%if close_master_ix>1
%    close_master_ix=[close_master_ix-1;close_master_ix];
%end

%dph_space=zeros(uw.n_edge,uw.n_ifg);
%dph_space_uw=zeros(uw.n_edge,uw.n_ifg);

%n_pad=2^ceil(log(n_image*2)/log(2));
%dph_padded=zeros(1,n_pad);
%smooth_win=gausswin(n_pad,alpha)';
%dph_smooth_uw=zeros(1,n_image-1);
%dph_noise=zeros(uw.n_edge,1);
dph_noise=zeros(ui.n_edge,uw.n_ifg,'single');
dph_space=((uw.ph(ui.edges(:,3),:).*conj(uw.ph(ui.edges(:,2),:))));
ph_noise=angle(uw.ph.*conj(uw.ph_lowpass));
clear uw
dph_noise_sf=((ph_noise(ui.edges(:,3),:)-(ph_noise(ui.edges(:,2),:))));
%dph_space_nc=angle(dph_space); 
%dph_space_sf=angle(dph_space.*exp(-j*dph_noise_sf)); % subtract noise estimate
%ix=std(dph_space_sf,0,2)<std(dph_space_nc,0,2);
%dph_space_nc(ix,:)=dph_space_sf(ix,:);
%dph_noise(ix,:)=dph_noise_sf(ix,:);
dph_space=dph_space./abs(dph_space);


if strcmpi(unwrap_method,'2D')
    save('uw_space_time','dph_space');
elseif strcmpi(unwrap_method,'3D_NO_DEF')
    dph_noise=angle(dph_space);
    dph_space_uw=angle(dph_space);
    save('uw_space_time','dph_space','dph_space_uw','dph_noise');
else
    %global G n dph x

    %G=zeros(n_ifg,n_image-1);
    %for i=1:n_ifg
    %    G(i,ifgday_ix(i,1):ifgday_ix(i,2)-1)=1;
    %end
    G=zeros(n_ifg,n_image);
    for i=1:n_ifg
        G(i,ifgday_ix(i,1))=-1;
        G(i,ifgday_ix(i,2))=1;
    end

    nzc_ix=sum(abs(G))~=0; % non-zero column index
    day=day(nzc_ix);
    G=G(:,nzc_ix);
    %n=n_image;
    n=size(G,2);
    
    %x=(0:n-1)'; % use sequence for smoothing
    x=(day-day(1))*(n-1)/(day(end)-day(1)); % use dates for smoothing
    
    %rand_ix=unique(round(rand(10000,1)*uw.n_edge))';
    %small_defo=pi/2;


    %small_ix=sum(abs(angle(dph_space))>small_defo,2)==0; % apparently signal is small
    %dph_space_series=[zeros(1,sum(small_ix));(G(:,2:n_image)\double(angle(dph_space(small_ix,:)))')];
    
    %dph_space_series=[zeros(1,uw.n_edge);double(W*G(:,2:n_image))\double(W*angle(dph_space)')];
    % No weighting needed here because model phase values include noise that's
    % common to all ifgs, and only the noise specific to each ifg is true noise
    % for this inversion, which is uncorrelated.  
    %dph_space_series=[zeros(1,ui.n_edge);double(G(:,2:n_image))\double(angle(dph_space))'];
    dph_space_series=[zeros(1,ui.n_edge);double(G(:,2:end))\double(angle(dph_space))'];
    dph_smooth_series=zeros(size(G,2),ui.n_edge,'single');

    for i1=1:n
        time_diff_sq=(day(i1)-day).^2;
        weight_factor=exp(-time_diff_sq/2/time_win^2);
        weight_factor=weight_factor/sum(weight_factor);
%        dph_smooth_series(i1,small_ix)=sum(dph_space_series.*repmat(weight_factor,1,sum(small_ix)));
        dph_smooth_series(i1,:)=sum(dph_space_series.*repmat(weight_factor,1,ui.n_edge));
    end
    %clear dph_space_series

    dph_smooth_ifg=(G*dph_smooth_series)';
    dph_noise=angle(dph_space.*exp(-j*dph_smooth_ifg));


    if strcmpi(unwrap_method,'3D_SMALL_DEF')|...
       strcmpi(unwrap_method,'3D_QUICK')
      %not_small_ix1=find(max(abs(angle(dph_space.*exp(-j*G*dph_space_series).')),[],2)>4)'; % inversion not trustworthy
      not_small_ix=find(std(dph_noise,0,2)>1.3)';
      %not_small_ix=unique([not_small_ix1,not_small_ix2]);
      fprintf('Ignoring %d edges. Elapsed time = %d s\n',length(not_small_ix),round(toc))
      dph_noise(not_small_ix,:)=nan;
    else
      m_minmax=repmat([-pi,pi],5,1).*repmat([0.5;0.25;1;0.25;1],1,2);
      anneal_opts=[1;15;0;0;0;0;0];
      %covm=cov(angle(dph_space)); % biased estimate of covariance
      covm=cov((dph_noise_sf)); % estimate of covariance
      [W,P]=chol(inv(covm)); % weighting matrix 
      if P~=0
          W=diag(1./sqrt(diag(covm)));
      end
      %max_noise=pi/2;
      %not_small_ix1=find(max(abs(angle(dph_space.*exp(-j*G*dph_space_series).')),[],2)>max_noise)'; % inversion not trustworthy
    %figure
    %hist(std(dph_noise,0,2),100)
      not_small_ix=find(std(dph_noise,0,2)>1)';
      %not_small_ix=unique([not_small_ix1,not_small_ix2]);
      fprintf('Performing more complex smoothing on %d edges. Elapsed time = %d s\n',length(not_small_ix),round(toc))

    %not_small_ix=find(~small_ix)';
      n_proc=0;
      for i=not_small_ix
        dph=angle(dph_space(i,:))';
        dph_smooth_series(:,i)=uw_sb_smooth_unwrap(m_minmax,anneal_opts,G,W,dph,x);
    %   dph_smooth_step=m(1)+m(2)*sin(2*pi/n*[0:n-1]'-m(3))+m(4)*sin(4*pi/n*[0:n-1]'-m(5));
        %dph_smooth_series(:,i)=m(1)*x+m(2)*n/2*sin(2*pi/n*x-m(3))+m(4)*n/2*sin(4*pi/n*x-m(5));

        n_proc=n_proc+1;
        if round(n_proc/1000)==n_proc/1000
            save('uw_unwrap_time','G','dph_space','dph_smooth_series');
            fprintf('%d edges of %d reprocessed. Elapsed time = %d s\n',n_proc,length(not_small_ix),round(toc))
        end
      end
      dph_smooth_ifg=(G*dph_smooth_series)';
      dph_noise=angle(dph_space.*exp(-j*dph_smooth_ifg));
    end
    dph_space_uw=dph_smooth_ifg+dph_noise;

    %[U,S,V]=svd(dph_space_uw(rand_ix,:));
    %SS=zeros(size(S));
    %SS(1:n_comp,1:n_comp)=S(1:n_comp,1_n_comp);
    %dph_smooth_ifg=U*SS*V';
    %dph_space_uw(rand_ix,:)=dph_smooth_ifg+angle(dph_space(i,:).*exp(-j*dph_smooth_ifg));

    %n_comp=5;
    %p_comp=S(1:n_comp,1:n_comp)*V(:,1:n_comp)';
    %m2_minmax=[min(U(:,1:n_comp))'*2,max(U(:,1:n_comp))'*2];

    %anneal_opts2=[1;10;0;0;0;0;1];

    %for i=1:uw.n_edge
    %    dph_space(i,:)=(uw.ph(uw.edges(i,3),:).*conj(uw.ph(uw.edges(i,2),:)));
    %    dph_space(i,:)=dph_space(i,:)./abs(dph_space(i,:));
    %    dph=dph_space(i,:);
    %    m2=anneal(@(x) sum(abs(angle(dph.*exp(-j*x'*p_comp)))),m2_minmax,anneal_opts2);
    %    dph_smooth_ifg=
    %    dph_noise(i,:)=angle(dph_space(i,:).*exp(-j*dph_smooth_ifg));
    %    dph_space_uw(i,:)=dph_smooth_ifg+angle(dph_space(i,:).*exp(-j*dph_smooth_ifg));


    %    if round(i/100)==i/100
    %        fprintf('%d edges of %d processed. Elapsed time = %d s\n',i,uw.n_edge,round(toc))
    %    end
    %end



    save('uw_space_time','dph_space','dph_space_uw','dph_noise','G','dph_smooth_series');
end
