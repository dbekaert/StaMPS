function []=PS_correct_phase
% PS_CORRECT_PHASE() correct phase from estimate of look angle error
%
%   Andy Hooper, June 2006
%
%   ==========================================================
%   07/2006 AH: Use specific bperp for correction
%   09/2006 AH: add small baselines 
%   ==========================================================

fprintf('Correcting phase for look angle error...\n')


small_baseline_flag=getparm('small_baseline_flag',1);

load psver
psname=['ps',num2str(psver)];
phname=['ph',num2str(psver)];
pmname=['pm',num2str(psver)];
rcname=['rc',num2str(psver)];
bpname=['bp',num2str(psver)];


ps=load(psname);
pm=load(pmname);
bp=load(bpname);


if exist([phname,'.mat'],'file')
    phin=load(phname);
    ph=phin.ph;
    clear phin
else
    ph=ps.ph;
end

K_ps=single(pm.K_ps);
C_ps=single(pm.C_ps);
master_ix=sum(ps.master_day>ps.day)+1;

if strcmpi(small_baseline_flag,'y')
    %ph_rc=ph.*exp(-j*(repmat(K_ps,1,ps.n_ifg).*bp.bperp_mat));  % subtract range error 
    ph_rc=ph.*exp(-j*(repmat(K_ps,1,ps.n_ifg).*bp.bperp_mat));  % subtract range error 
    save(rcname,'ph_rc'); 
else    
    bperp_mat=[bp.bperp_mat(:,1:ps.master_ix-1),zeros(ps.n_ps,1,'single'),bp.bperp_mat(:,ps.master_ix:end)];
    ph_rc=ph.*exp(-j*(repmat(K_ps,1,ps.n_ifg).*bperp_mat+repmat(C_ps,1,ps.n_ifg)));  % subtract range error and master noise
    ph_reref=[single(pm.ph_patch(:,1:master_ix-1)),ones(ps.n_ps,1,'single'),single(pm.ph_patch(:,master_ix:end))];
    save(rcname,'ph_rc','ph_reref'); 
end

    
    

