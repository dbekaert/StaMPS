function []=PS_calc_ifg_std
% PS_CALC_IFG_STD() calculate std for each ifg
%
%   Andy Hooper, June 2006
%
%   ======================================================
%   09/2006 AH: small baselines added 
%   ======================================================

fprintf('Estimating noise standard deviation...\n')

small_baseline_flag=getparm('small_baseline_flag');

load psver
psname=['ps',num2str(psver)];
phname=['ph',num2str(psver)];
pmname=['pm',num2str(psver)];
bpname=['bp',num2str(psver)];
ifgstdname=['ifgstd',num2str(psver)];

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

n_ps=length(ps.xy);
master_ix=sum(ps.master_day>ps.day)+1;

if strcmpi(small_baseline_flag,'y')
    ph_diff=angle(ph.*conj(pm.ph_patch).*exp(-j*(repmat(pm.K_ps,1,ps.n_ifg).*bp.bperp_mat)));    
else
    bperp_mat=[bp.bperp_mat(:,1:ps.master_ix-1),zeros(ps.n_ps,1,'single'),bp.bperp_mat(:,ps.master_ix:end)];
    ph_patch=[pm.ph_patch(:,1:master_ix-1),ones(n_ps,1),pm.ph_patch(:,master_ix:end)];
    ph_diff=angle(ph.*conj(ph_patch).*exp(-j*(repmat(pm.K_ps,1,ps.n_ifg).*bperp_mat+repmat(pm.C_ps,1,ps.n_ifg))));
end

ifg_std=[sqrt(sum(ph_diff.^2)/n_ps)*180/pi]'


save(ifgstdname,'ifg_std'); 
    
    

