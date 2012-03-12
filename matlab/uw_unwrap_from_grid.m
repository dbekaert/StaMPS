function [ph_uw,msd]=uw_unwrap_from_grid(ph,xy,pix_size)
%UW_UNWRAP_FROM_GRID unwrap PS from unwrapped gridded ifgs 
%
%   Andy Hooper, June 2007
%
% ============================================================
% 03/2012 AH: Allow for all zero wrapped phase values
% 03/2012 AH: Allow for non-complex wrapped phase
% ============================================================

fprintf('Unwrapping from grid...\n')

uw=load('uw_grid');
uu=load('uw_phaseuw');

[n_ps,n_ifg]=size(ph);
gridix=zeros(size(uw.nzix));
gridix(uw.nzix)=[1:uw.n_ps];
    
ph_uw=zeros(n_ps,n_ifg,'single');

for i=1:n_ps
    ix=gridix(uw.grid_ij(i,1),uw.grid_ij(i,2));
    if ix==0
        ph_uw(i,:)=nan; % wrapped phase values were zero
    else
        ph_uw_pix=uu.ph_uw(ix,:);
        if isreal(ph)
          ph_uw(i,:)=ph_uw_pix+angle(exp(1i*(ph(i,:)-ph_uw_pix)));
        else
          ph_uw(i,:)=ph_uw_pix+angle(ph(i,:).*exp(-1i*ph_uw_pix));
        end
    end
end

msd=uu.msd;
