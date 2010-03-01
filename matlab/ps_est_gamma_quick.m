function []=PS_est_gamma_quick(restart_flag)
% PS_EST_GAMMA_QUICK estimate coherence of PS cands
%   PS_EST_GAMMA_QUICK(RESTART_FLAG) set R
%   RESTART_FLAG=1 restarts from previous values
%   RESTART_FLAG=2 restarts from previous values, but only calculates patch values
%   
%
%   Andy Hooper, June 2006
%
%   ==========================================================
%   09/2006 AH: short baseline added,  
%   09/2006 AH: unwrapped phase loaded from separate workspace  
%   09/2006 AH: restart processing fixed 
%   10/2006 AH: convergence criteria changed
%   04/2007 AH: number of wraps no longer rounded
%   12/2008 AH: avoid divide by zero for zero phase values
%   ==========================================================

fprintf('Estimating gamma for candidate pixels...\n')


if nargin<1
   restart_flag=0;
end

rho = 830000; % mean range - need only be approximately correct
n_rand=300000; % number of simulated random phase pixels

grid_size=getparm('filter_grid_size',1);
filter_weighting=getparm('filter_weighting',1);
n_win=getparm('clap_win',1);
low_pass_wavelength=getparm('clap_low_pass_wavelength',1);
clap_alpha=getparm('clap_alpha',1);
clap_beta=getparm('clap_beta',1);
max_topo_err=getparm('max_topo_err',1);
lambda=getparm('lambda',1);
gamma_change_convergence=getparm('gamma_change_convergence',1);
small_baseline_flag=getparm('small_baseline_flag',1);

if strcmpi(small_baseline_flag,'y')
    low_coh_thresh=15; % equivalent to coh of 15/100
else
    low_coh_thresh=31; % equivalent to coh of 31/100
end

freq0=1/low_pass_wavelength;
freq_i=-(n_win)/grid_size/n_win/2:1/grid_size/n_win:(n_win-2)/grid_size/n_win/2;
butter_i=1./(1+(freq_i/freq0).^(2*5));
low_pass=butter_i'*butter_i;
low_pass=fftshift(low_pass);

load psver
psname=['ps',num2str(psver)];
phname=['ph',num2str(psver)];
bpname=['bp',num2str(psver)];
laname=['la',num2str(psver),'.mat'];
pmname=['pm',num2str(psver),'.mat'];
daname=['da',num2str(psver),'.mat'];

ps=load(psname);

bp=load(bpname);

if exist(daname,'file')
    da=load(daname);
    D_A=da.D_A;
    clear da
else
    D_A=ones(ps.n_ps,1);
end

if exist([phname,'.mat'],'file')
    phin=load(phname);
    ph=phin.ph;
    clear phin
else
    ph=ps.ph;
end

[null_i,null_j]=find(ph==0);
null_i=unique(null_i);
good_ix=logical(ones(ps.n_ps,1));
good_ix(null_i)=0;

if strcmpi(small_baseline_flag,'y')
    bperp=ps.bperp;
    n_ifg=ps.n_ifg;
    n_image=ps.n_image;
    n_ps=ps.n_ps;
    ifgday_ix=ps.ifgday_ix;
    xy=ps.xy;
else
    ph=ph(:,[1:ps.master_ix-1,ps.master_ix+1:end]);
    bperp=ps.bperp([1:ps.master_ix-1,ps.master_ix+1:end]);
    n_ifg=ps.n_ifg-1;
    n_ps=ps.n_ps;
    xy=ps.xy;
end
clear ps

A=abs(ph);
A=single(A);
A(A==0)=1; % avoid divide by zero
ph=ph./A;

if exist(laname,'file')
    la=load(laname);
    inc_mean=mean(la.la)+0.052; % incidence angle approx equals look angle + 3 deg
    clear la
else
    inc_mean=21*pi/180 % guess the incidence angle
end

max_K=max_topo_err/(lambda*rho*sin(inc_mean)/4/pi);
bperp_range=max(bperp)-min(bperp);
n_trial_wraps=(bperp_range*max_K/(2*pi))


if restart_flag > 0
    %disp(['Restarting: iteration #',num2str(i_loop),' step_number=',num2str(step_number)])
    fprintf('Restarting previous run...\n')
    load(pmname)
    weighting_save=weighting;
    if ~exist('gamma_change_save','var')
        gamma_change_save=1;
    end


else
    fprintf('Initialising random distribution...\n')
    rand('state',2005)
    % determine distribution for random phase

    if strcmpi(small_baseline_flag,'y')
         rand_image=2*pi*rand(n_rand,n_image);
         rand_ifg=zeros(n_rand,n_ifg);
         for i=1:n_ifg
             rand_ifg(:,i)=rand_image(:,ifgday_ix(i,2))-rand_image(:,ifgday_ix(i,1));
         end
         clear rand_image
    else
         rand_ifg=2*pi*rand(n_rand,n_ifg);
    end
    for i=n_rand:-1:1      
        [K_r,C_r,coh_r]=ps_topofit(exp(j*rand_ifg(i,:)),bperp,n_trial_wraps,'n');
        coh_rand(i)=coh_r(1);
    end
    clear rand_ifg
    coh_bins=[0.005:0.01:0.995];
    Nr=hist(coh_rand,coh_bins); % distribution of random phase points
    i=length(Nr);
    while Nr(i)==0
        i=i-1;
    end
    Nr_max_nz_ix=i;

    step_number=1;
    K_ps=zeros(n_ps,1);
    C_ps=zeros(n_ps,1);
    coh_ps=zeros(n_ps,1);
    coh_ps_save=zeros(n_ps,1);
    N_opt=zeros(n_ps,1);
    ph_res=zeros(n_ps,n_ifg,'single');
    ph_patch=zeros(size(ph),'single');
    N_patch=zeros(n_ps,1);
    grid_ij(:,1)=ceil((xy(:,3)-min(xy(:,3))+1e-6)/grid_size);
    grid_ij(grid_ij(:,1)==max(grid_ij(:,1)),1)=max(grid_ij(:,1))-1;
    grid_ij(:,2)=ceil((xy(:,2)-min(xy(:,2))+1e-6)/grid_size);
    grid_ij(grid_ij(:,2)==max(grid_ij(:,2)),2)=max(grid_ij(:,2))-1;
    i_loop=1;
    weighting=1./D_A; 
    weighting_save=weighting;
    gamma_change_save=0;

end

n_i=max(grid_ij(:,1));
n_j=max(grid_ij(:,2));


fprintf('%d PS candidates to process\n',n_ps)
xy(:,1)=[1:n_ps]'; % assumption that already sorted in ascending column 3 (y-axis) order
loop_end_sw=0;
n_high_save=0;

while loop_end_sw==0
  %if step_number==1     % check in case restarting and step 1 already completed
    fprintf('\niteration #%d\n',i_loop)
    fprintf('Calculating patch phases...\n')

    ph_grid=zeros(n_i,n_j,n_ifg,'single');
    ph_filt=ph_grid;
    
    for i=1:n_ps
        ph_grid(grid_ij(i,1),grid_ij(i,2),:)=ph_grid(grid_ij(i,1),grid_ij(i,2),:)+shiftdim(ph(i,:),-1)*weighting(i);
    end
    
    for i=1:n_ifg
        ph_filt(:,:,i)=clap_filt(ph_grid(:,:,i),clap_alpha,clap_beta,n_win*0.75,n_win*0.25,low_pass);
    end
        
    for i=1:n_ps
        ph_patch(i,1:n_ifg)=squeeze(ph_filt(grid_ij(i,1),grid_ij(i,2),:));
    end
        
    clear ph_filt
    ix=ph_patch~=0;
    ph_patch(ix)=ph_patch(ix)./abs(ph_patch(ix));
    %ph_patch=ph_patch./abs(ph_patch);
  %end % end-if step_number
  
%%%%%%%%%%%%% Now estimate topo error %%%%%%%%%%%%%%%%%%%%%
    if restart_flag<2
    
        fprintf('Estimating topo error...\n')
        step_number=2;

        for i=1:n_ps
            psdph=ph(i,:).*conj(ph_patch(i,:));
            if sum(psdph==0)==0  % insist on a non-null value in every ifg
                [Kopt,Copt,cohopt,ph_residual]=ps_topofit(psdph,bp.bperp_mat(i,:)',n_trial_wraps,'n');
                K_ps(i)=Kopt(1);
                C_ps(i)=Copt(1);
                coh_ps(i)=cohopt(1);
                N_opt(i)=length(Kopt);
                ph_res(i,:)=angle(ph_residual);
            else
                K_ps(i)=nan;
                coh_ps(i)=0;
            end
            if i/100000==floor(i/100000)
                fprintf('%d PS processed\n',i)
            end
        end
        
        
        if strcmpi(filter_weighting,'P-square')
            Na=hist(coh_ps,coh_bins);
            Nr=Nr*sum(Na(1:low_coh_thresh))/sum(Nr(1:low_coh_thresh)); % scale random distribution to actual, using low coh values
            Na(Na==0)=1; % avoid divide by zero
            Prand=Nr./Na;
            Prand(1:low_coh_thresh)=1;
            Prand(Nr_max_nz_ix+1:end)=0;
            Prand(Prand>1)=1;
            Prand=filter(gausswin(7),1,[ones(1,7),Prand])/sum(gausswin(7));
            Prand=Prand(8:end);
            Prand=interp([1,Prand],10); % interpolate to 100 samples
            Prand=Prand(1:end-9);
            Prand_ps=Prand(round(coh_ps*1000)+1)';
            weighting=(1-Prand_ps).^2;
        else
            %ph_n=angle(ph_res.*repmat(conj(sum(ph_res,2)),1,n_ifg)); % subtract mean, take angle
            %sigma_n=std(A.*sin(ph_n),0,2); % noise
            
            g=mean(A.*cos(ph_res),2); % signal
            sigma_n=sqrt(0.5*(mean(A.^2,2)-g.^2));
            %snr=(g./sigma_n).^2;

            weighting(sigma_n==0)=0;
            weighting(sigma_n~=0)=g(sigma_n~=0)./sigma_n(sigma_n~=0); % snr
        end

        weighting_rms=sqrt(sum((weighting-weighting_save).^2)/n_ps);
        weighting_save=weighting;


        step_number=1;

        %if i_loop==1
        %    figure
        %    subplot(2,1,1)
        %    hist(weighting,100)
        %    subplot(2,1,2)
        %    hist(coh_ps,100)
        %end

        gamma_change_rms=sqrt(sum((coh_ps-coh_ps_save).^2)/n_ps);
        gamma_change_change=gamma_change_rms-gamma_change_save
        gamma_change_save=gamma_change_rms;
        coh_ps_save=coh_ps;

        gamma_change_convergence=getparm('gamma_change_convergence',1);
        if gamma_change_change<0&abs(gamma_change_change)<gamma_change_convergence
            %figure
            %subplot(2,1,1)
            %hist(weighting,100)
            %subplot(2,1,2)
            %hist(coh_ps,100)
            loop_end_sw=1;
        else
            i_loop=i_loop+1;
        end
    else
        loop_end_sw=1;
    end
    
save(pmname,'ph_patch','K_ps','C_ps','coh_ps','N_opt','ph_res','step_number','ph_grid','n_trial_wraps','grid_ij','grid_size','low_pass','i_loop','weighting','Nr','Nr_max_nz_ix','coh_bins','coh_ps_save','gamma_change_save') 
end

