function [ph_out]=clap_filt_patch(ph,alpha,beta,low_pass)
%CLAP_FILT_PATCH Combined Low-pass Adaptive Phase filtering on 1 patch
%   [ph_out]=clap_filt_patch(ph,alpha,beta)
%
%
%   Andy Hooper, June 2006



if nargin<2
    alpha=0.5;
end

if nargin<3
    beta=0.1;
end

if nargin<4 | isempty(low_pass)
    low_pass=zeros(size(ph));
end


ph(isnan(ph))=0;
B=gausswin(7)*gausswin(7)';

ph_fft=fft2(ph);
H=abs(ph_fft);
H=ifftshift(filter2(B,fftshift(H))); % smooth response
meanH=median(H(:));
if meanH~=0
    H=H/meanH;
end
H=H.^alpha;
H=H-1; % set all values under median to zero
H(H<0)=0; % set all values under median to zero
G=H*beta+low_pass;
ph_out=ifft2(ph_fft.*G);
