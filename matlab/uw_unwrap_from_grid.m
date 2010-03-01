function [ph_uw,msd]=uw_unwrap_from_grid(ph,xy,pix_size)
%UW_UNWRAP_FROM_GRID unwrap PS from unwrapped gridded ifgs 
%
%   Andy Hooper, June 2007
%
fprintf('Unwrapping from grid...\n')

uw=load('uw_grid');
uu=load('uw_phaseuw');

[n_ps,n_ifg]=size(ph);
gridix=zeros(size(uw.nzix));
gridix(uw.nzix)=[1:uw.n_ps];
    
ph_uw=zeros(n_ps,n_ifg);

for i=1:n_ps
    ph_uw_pix=uu.ph_uw(gridix(uw.grid_ij(i,1),uw.grid_ij(i,2)),:);
    ph_uw(i,:)=ph_uw_pix+angle(ph(i,:).*exp(-j*ph_uw_pix));
end

msd=uu.msd;
