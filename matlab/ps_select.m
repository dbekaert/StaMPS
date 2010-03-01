function []=ps_select(reest_flag,plot_flag)
% PS_SELECT select PS based on gamma and D_A
%   PS_SELECT(REEST_FLAG,PLOT_FLAG)  
%   REEST_FLAG = 1 skip reestimation completely (default 0)
%   REEST_FLAG = 2 reuse previously reestimated gamma (default 0)
%   REEST_FLAG = 3 reestimate gamma for all candidate pixels (default 0)
%   PLOT_FLAG = 1 for plots (default 0)
%
%   Andy Hooper, June 2006
%
%   ======================================================================
%   07/2006 AH: Use specific bperp for re-estimation of gamma
%   08/2006 AH: Load workspaces into structures
%   09/2006 AH: Check if no gamma threshold meets criteria
%   09/2006 AH: add small baselines
%   09/2009 AH: add option to select on density instead of percentage
%   01/2010 KS: fixed badly conditioned polyfit errors by scaling and
%               centering
%   02/2010 AH: Only bin by D_A if enough candidate pixels
%   ======================================================================

fprintf('Selecting stable-phase pixels...\n')

if nargin<1
    reest_flag=0;  
end

if nargin<2
    plot_flag=0;  
end

psver=getpsver;
if psver>1
   setpsver(1)
end

clap_alpha=getparm('clap_alpha',1);
clap_beta=getparm('clap_beta',1);
n_win=getparm('clap_win',1);
select_method=getparm('select_method',1);
if strcmpi(select_method,'PERCENT')
    max_percent_rand=getparm('percent_rand',1);
else
    max_density_rand=getparm('density_rand',1);
end
gamma_stdev_reject=getparm('gamma_stdev_reject',1);
small_baseline_flag=getparm('small_baseline_flag',1);


if strcmpi(small_baseline_flag,'y')
    low_coh_thresh=15; % equivalent to coh of 15/100
else
    low_coh_thresh=31; % equivalent to coh of 31/100
end

load psver

psname=['ps',num2str(psver)];
phname=['ph',num2str(psver)];
pmname=['pm',num2str(psver)];
selectname=['select',num2str(psver)];
daname=['da',num2str(psver)];
bpname=['bp',num2str(psver)];

ps=load(psname);

if exist([phname,'.mat'],'file')
    phin=load(phname);
    ph=phin.ph;
    clear phin
else
    ph=ps.ph;
end

if strcmpi(small_baseline_flag,'y')
    bperp=ps.bperp;
    n_ifg=ps.n_ifg;
    n_ps=ps.n_ps;
    xy=ps.xy;
else
    master_ix=sum(ps.master_day>ps.day)+1;
    ph=ph(:,[1:ps.master_ix-1,ps.master_ix+1:end]);
    bperp=ps.bperp([1:ps.master_ix-1,ps.master_ix+1:end]);
    n_ifg=ps.n_ifg-1;
    n_ps=ps.n_ps;
    xy=ps.xy;
end
clear ps

pm=load(pmname);
if exist([daname,'.mat'],'file')
    da=load(daname);
    D_A=da.D_A;
    clear da
else
    D_A=[];
end

if ~isempty(D_A) & size(D_A,1)>=10000
    % chunk up PSC
    D_A_sort=sort(D_A);
    if size(D_A,1)>=50000
      bin_size=10000;
    else
      bin_size=2000;
    end
    D_A_max=[0;D_A_sort(bin_size:bin_size:end-bin_size);D_A_sort(end)];

else
    D_A_max=[0;1];
    D_A=ones(size(pm.coh_ps));
end

if ~strcmpi(select_method,'PERCENT')
    patch_area=prod(max(xy(:,2:3),[],1)-min(xy(:,2:3),[],1))/1e6; % in km
    max_percent_rand=max_density_rand*patch_area/(length(D_A_max)-1);
end

min_coh=zeros(length(D_A_max)-1,1);
D_A_mean=zeros(size(D_A_max,1)-1,1);
Nr_dist=pm.Nr;

if reest_flag==3 % reestimate for all candidate pixels
    coh_thresh=0;
    coh_thresh_coeffs=[];
else
  for i=1:length(D_A_max)-1
    coh_chunk=pm.coh_ps(D_A>D_A_max(i) & D_A<=D_A_max(i+1));
    D_A_mean(i)=mean(D_A(D_A>D_A_max(i) & D_A<=D_A_max(i+1)));
    coh_chunk=coh_chunk(coh_chunk~=0);  % discard PSC for which coherence was not calculated
    Na=hist(coh_chunk,pm.coh_bins);
    Nr=Nr_dist*sum(Na(1:low_coh_thresh))/sum(Nr_dist(1:low_coh_thresh));
    if i==length(D_A_max)-1 & plot_flag==1 
        figure
        plot(pm.coh_bins,Na,'g')
        hold on
        plot(pm.coh_bins,Nr,'r')
        legend('data','random')
        title('Before Gamma Reestimation')
    end
    Na(Na==0)=1; % avoid divide by zero
    if strcmpi(select_method,'PERCENT')
        percent_rand=fliplr(cumsum(fliplr(Nr))./cumsum(fliplr(Na))*100);
    else
        percent_rand=fliplr(cumsum(fliplr(Nr))); % absolute number
    end
    ok_ix=find(percent_rand<max_percent_rand);
    if isempty(ok_ix)
        min_coh(i)=1; % no threshold meets criteria
    else
        min_fit_ix=min(ok_ix)-3;
        if min_fit_ix<=0
            %min_coh(i)=low_coh_thresh/100;
            min_coh(i)=nan;
        else
            max_fit_ix=min(ok_ix)+2;
            max_fit_ix(max_fit_ix>100)=100; % make sure not higher than length of percent_rand
            [p,S,mu]=polyfit(percent_rand(min_fit_ix:max_fit_ix),[min_fit_ix*0.01:0.01:max_fit_ix*0.01],3); % KS
            min_coh(i)=polyval(p,max_percent_rand,[],mu); % KS
%             p=polyfit(percent_rand(min_fit_ix:max_fit_ix),[min_fit_ix*0.01:0.01:max_fit_ix*0.01],3);
%             min_coh(i)=polyval(p,max_percent_rand);
        end
    end
  end
  nonnanix=~isnan(min_coh);
  if sum(nonnanix)<1
    error('could not set threshold - try setting percent_rand lower')
  end
  min_coh=min_coh(nonnanix);
  D_A_mean=D_A_mean(nonnanix);
  if size(min_coh,1)>1
    coh_thresh_coeffs=polyfit(D_A_mean,min_coh,1);  % fit polynomial to the curve
    if coh_thresh_coeffs(1)>0 % positive slope (as expected)
        coh_thresh=polyval(coh_thresh_coeffs,D_A);
    else % unable to ascertain correct slope
        coh_thresh=polyval(coh_thresh_coeffs,0.35); % set an average threshold for all D_A
        coh_thresh_coeffs=[];
    end
  else
    coh_thresh=min_coh;
    coh_thresh_coeffs=[];
  end
end

fprintf('Initial gamma threshold: %.3f at D_A=%.2f to %.3f at D_A=%.2f\n',min(coh_thresh),min(D_A),max(coh_thresh),max(D_A))

if  plot_flag==1 
    figure
    plot(D_A_mean,min_coh,'*')
    hold on
    if ~isempty(coh_thresh_coeffs)
        plot(D_A_mean,polyval(coh_thresh_coeffs,D_A_mean),'r')
    end
    ylabel('\gamma_{thresh}')
    xlabel('D_A')
end

ix=find(pm.coh_ps>coh_thresh); % select those below threshold
n_ps=length(ix);
fprintf('%d PS selected initially\n',n_ps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% reject part-time PS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gamma_stdev_reject>0
    ph_res_cpx=exp(j*pm.ph_res);
    coh_std=zeros(length(ix),1);
    for i=1:length(ix)
        coh_std(i)=std(bootstrp(100,@(ph) abs(sum(ph))/length(ph), ph_res_cpx(ix(i),:)));
    end
    clear ph_res_cpx
    ix=ix(coh_std<gamma_stdev_reject);
    n_ps=length(ix);
    fprintf('%d PS left after pps rejection\n',n_ps)
end


if reest_flag~=1  

  if reest_flag~=2 % reestimate coh with the PS removed from filtered patch
    pm=rmfield(pm,{'ph_res'});
    pm=rmfield(pm,{'ph_patch'});
    ph_patch2=zeros(n_ps,n_ifg,'single');
    ph_res2=zeros(n_ps,n_ifg,'single');
    ph=ph(ix,:);


    if length(coh_thresh)>1
        coh_thresh=coh_thresh(ix);
    end

    n_i=max(pm.grid_ij(:,1));
    n_j=max(pm.grid_ij(:,2));
    K_ps2=zeros(n_ps,1);
    C_ps2=zeros(n_ps,1);
    coh_ps2=zeros(n_ps,1);
    ph_filt=zeros(n_win,n_win,n_ifg);

    for i=1:n_ps
        ps_ij=pm.grid_ij(ix(i),:);
        i_min=max(ps_ij(1)-n_win/2,1);
        i_max=i_min+n_win-1;
        if i_max>n_i
            i_min=i_min-i_max+n_i;
            i_max=n_i;
        end
        j_min=max(ps_ij(2)-n_win/2,1);
        j_max=j_min+n_win-1;
        if j_max>n_j
            j_min=j_min-j_max+n_j;
            j_max=n_j;
        end
        ps_bit_i=ps_ij(1)-i_min+1;
        ps_bit_j=ps_ij(2)-j_min+1;

        ph_bit=pm.ph_grid(i_min:i_max,j_min:j_max,:);

        ph_bit(ps_bit_i,ps_bit_j,:)=0;

        for i_ifg=1:n_ifg
            ph_filt(:,:,i_ifg)=clap_filt_patch(ph_bit(:,:,i_ifg),clap_alpha,clap_beta,pm.low_pass);
        end

        ph_patch2(i,:)=squeeze(ph_filt(ps_bit_i,ps_bit_j,:));
        if i/10000==floor(i/10000)
            fprintf('%d patches re-estimated\n',i)
        end

    end
    
    pm=rmfield(pm,{'ph_grid'});
    bp=load(bpname);
    bperp_mat=bp.bperp_mat(ix,:);
    clear bp

    for i=1:n_ps
        psdph=ph(i,:).*conj(ph_patch2(i,:));
        if sum(psdph==0)==0  % insist on a non-null value in every ifg
            psdph=psdph./abs(psdph);
            [Kopt,Copt,cohopt,ph_residual]=ps_topofit(psdph,bperp_mat(i,:)',pm.n_trial_wraps,'n');
            K_ps2(i)=Kopt(1);
            C_ps2(i)=Copt(1);
            coh_ps2(i)=cohopt(1);
            ph_res2(i,:)=angle(ph_residual);
        else
            K_ps2(i)=nan;
            coh_ps2(i)=nan;
        end
        if i/10000==floor(i/10000)
            fprintf('%d coherences re-estimated\n',i)
        end
    end
  else % reest_flag==2, use previously recalculated coh 
    sl=load(selectname);
    ix=sl.ix;
    coh_ps2=sl.coh_ps2;
    K_ps2=sl.K_ps2;
    C_ps2=sl.C_ps2;
    ph_res2=sl.ph_res2;
    ph_patch2=sl.ph_patch2;
  end 
    pm.coh_ps(ix)=coh_ps2;

    % calculate threshold again based on recalculated coh
    for i=1:length(D_A_max)-1
        coh_chunk=pm.coh_ps(D_A>D_A_max(i) & D_A<=D_A_max(i+1));
        D_A_mean(i)=mean(D_A(D_A>D_A_max(i) & D_A<=D_A_max(i+1)));
        coh_chunk=coh_chunk(coh_chunk~=0);  % discard PSC for which coherence was not calculated
        Na=hist(coh_chunk,pm.coh_bins);
        Nr=Nr_dist*sum(Na(1:low_coh_thresh))/sum(Nr_dist(1:low_coh_thresh));
        if i==length(D_A_max)-1 & plot_flag==1 
            figure
            plot(pm.coh_bins,Na,'g')
            hold on
            plot(pm.coh_bins,Nr,'r')
            legend('data','random')
            title('After Gamma Reestimation')
        end
        Na(Na==0)=1; % avoid divide by zero
        if strcmpi(select_method,'PERCENT')
            percent_rand=fliplr(cumsum(fliplr(Nr))./cumsum(fliplr(Na))*100);
        else
            percent_rand=fliplr(cumsum(fliplr(Nr)));
        end
        ok_ix=find(percent_rand<max_percent_rand);
        if isempty(ok_ix)
            min_coh(i)=1;
        else
            min_fit_ix=min(ok_ix)-3;
            if min_fit_ix<=0
                %min_coh(i)=low_coh_thresh/100;
                min_coh(i)=nan;
            else
                max_fit_ix=min(ok_ix)+2;
                max_fit_ix(max_fit_ix>100)=100; % make sure not higher than length of percent_rand
                [p,S,mu]=polyfit(percent_rand(min_fit_ix:max_fit_ix),[min_fit_ix*0.01:0.01:max_fit_ix*0.01],3); % KS
                min_coh(i)=polyval(p,max_percent_rand,[],mu); % KS
%                 p=polyfit(percent_rand(min_fit_ix:max_fit_ix),[min_fit_ix*0.01:0.01:max_fit_ix*0.01],3);
%                 min_coh(i)=polyval(p,max_percent_rand);
            end
        end
    end

    nonnanix=~isnan(min_coh);
    min_coh=min_coh(nonnanix);
    D_A_mean=D_A_mean(nonnanix);
    if size(min_coh,1)>1
        coh_thresh_coeffs=polyfit(D_A_mean,min_coh,1);  % fit polynomial to the curve
        if coh_thresh_coeffs(1)>0 % positive slope (as expected)
            coh_thresh=polyval(coh_thresh_coeffs,D_A(ix));
        else % unable to ascertain correct slope
            coh_thresh=polyval(coh_thresh_coeffs,0.35); % set an average threshold for all D_A
            coh_thresh_coeffs=[];
        end
    else
        coh_thresh=min_coh;
        coh_thresh_coeffs=[];
    end

    fprintf('Reestimation gamma threshold: %.3f at D_A=%.2f to %.3f at D_A=%.2f\n',min(coh_thresh),min(D_A),max(coh_thresh),max(D_A))


    bperp_range=max(bperp)-min(bperp);
    keep_ix=coh_ps2>coh_thresh & abs(pm.K_ps(ix)-K_ps2)<2*pi/bperp_range;
    clear pm
    fprintf('%d ps selected after re-estimation of coherence\n',sum(keep_ix))
    clear pm

else % reest_flag==1, skip re-estimation
    pm=rmfield(pm,{'ph_grid'});
    ph_patch2=pm.ph_patch(ix,:);
    ph_res2=pm.ph_res(ix,:);
    K_ps2=pm.K_ps(ix);
    C_ps2=pm.C_ps(ix);
    coh_ps2=pm.coh_ps(ix);
    keep_ix=true(size(ix));
end % end-if reest_flag

if  plot_flag==1 
    figure
    plot(D_A_mean,min_coh,'*')
    hold on
    if ~isempty(coh_thresh_coeffs)
        plot(D_A_mean,polyval(coh_thresh_coeffs,D_A_mean),'r')
    end
    ylabel('\gamma_{thresh}')
    xlabel('D_A')
end


save(selectname,'ix','keep_ix','ph_patch2','ph_res2','K_ps2','C_ps2','coh_ps2','coh_thresh','coh_thresh_coeffs','clap_alpha','clap_beta','n_win','max_percent_rand','gamma_stdev_reject','small_baseline_flag');
    

