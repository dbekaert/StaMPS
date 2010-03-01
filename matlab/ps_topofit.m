function [K0,C0,coh0,phase_residual]=ps_topofit(cpxphase,bperp,n_trial_wraps,plotflag)
%PS_TOPOFIT find best-fitting range error 
%   PS_TOPOFIT(cpxphase,bperp,n_trial_wraps,plotflag)
%
%   Andy Hooper, June 2006
%
%   ==========================================================
%   04/2007 AH: Added 64-bit machine compatibility
%   04/2007 AH: Tightened up max topo error processing
%   ==========================================================



if size(cpxphase,2)>1
   cpxphase=cpxphase.';
end

ix=cpxphase~=0;  % if signal of one image is 0, dph set to 0
cpxphase=cpxphase(ix);
bperp=bperp(ix);
n_ix=length(ix);
bperp_range=max(bperp)-min(bperp);

wphase=angle(cpxphase);

trial_mult=[-ceil(8*n_trial_wraps):ceil(8*n_trial_wraps)];
n_trials=length(trial_mult);
trial_phase=bperp/bperp_range*pi/4;
trial_phase_mat=exp(-j*trial_phase*trial_mult);
cpxphase_mat=repmat(cpxphase,1,n_trials);
phaser=trial_phase_mat.*cpxphase_mat;
phaser_sum=sum(phaser);
C_trial=angle(phaser_sum);
coh_trial=abs(phaser_sum)/sum(abs(cpxphase));

coh_diff=diff(coh_trial);
coh_max_ix=[];
for i=2:length(coh_diff)
    if coh_diff(i)<0 & coh_diff(i-1)>0
        coh_max_ix=[coh_max_ix,i];
    end
end 

coh_max=coh_trial(coh_max_ix); % maximum value of coherence

[dummy,coh_high_max_ix]=max(coh_trial); % only select highest

K0=pi/4/bperp_range*trial_mult(coh_high_max_ix);
C0=C_trial(coh_high_max_ix);
coh0=coh_trial(coh_high_max_ix);


% linearise and solve
resphase=cpxphase.*exp(-j*(K0*bperp)); % subtract approximate fit
offset_phase=sum(resphase);
resphase=angle(resphase*conj(offset_phase)); % subtract offset, take angle (unweighted)
weighting=abs(cpxphase); 
mopt=double(weighting.*bperp)\double(weighting.*resphase);
K0=K0+mopt;
phase_residual=cpxphase.*exp(-j*(K0*bperp)); 
mean_phase_residual=sum(phase_residual); 
C0=angle(mean_phase_residual);    % static offset (due to noise of master + average noise of rest)
coh0=abs(mean_phase_residual)/sum(abs(phase_residual)); 

if plotflag=='y'
    subplot(2,1,2)
    bvec=linspace(min(bperp),max(bperp),200);
	wphase_hat=angle(exp(j*(K0(1)*bvec+C0(1))));
	p=plot(bvec,(wphase_hat),'r');
	hold on
    set(p,'linewidth',2)
    p=plot(bperp,wphase,'bo');
    set(p,'linewidth',2)
	hold off
    set(gca,'ylim',[-pi,pi])
    set(gca,'fontsize',12,'fontweight','bold')
    ylabel('Wrapped Phase')
    xlabel('B_{\perp} (m)')
    subplot(2,1,1)
    plot(pi/4/bperp_range/4/pi*0.05656*trial_mult,coh_trial,'g')
    ylabel('\gamma_x')
    xlabel('\Delta \theta^{nc}_x (radians)')
    set(gca,'fontsize',12,'fontweight','bold')
    axis tight
    set(gca,'ylim',[0,1])

end
