function []=uw_unwrap_time(alpha,master_day)
%UW_UNWRAP_TIME smooth and unwrap phase diffs between neighboring data points in time
%
%   Andy Hooper, June 2006

if nargin<1
    help uw_unwrap_time
    error('not enough arguments')
end

if nargin<2
    master_day=0;
end

disp('Unwrapping in time...')

uw=load('uw_phase');

day_pos_ix=find(uw.day>0);
[Y,I]=min(uw.day(day_pos_ix));
close_master_ix=day_pos_ix(I);
if close_master_ix>1
    close_master_ix=[close_master_ix-1;close_master_ix];
end

dph_space=zeros(uw.n_edge,uw.n_ifg);
dph_space_uw=zeros(uw.n_edge,uw.n_ifg);

n_pad=2^ceil(log(uw.n_ifg*2)/log(2));
dph_padded=zeros(1,n_pad);
smooth_win=gausswin(n_pad,alpha)';
dph_smooth_uw=zeros(1,uw.n_ifg);
dph_noise=zeros(uw.n_edge,uw.n_ifg);

for i=1:uw.n_edge
    
    dph_space(i,:)=(uw.ph(uw.edgs(i,3),:).*conj(uw.ph(uw.edgs(i,2),:)));
    dph_padded(1:uw.n_ifg)=dph_space(i,:);
    dph_fft=fftshift(fft(dph_padded));
    dph_smooth=ifft(ifftshift(dph_fft.*smooth_win));
    dph_noise(i,:)=angle(dph_space(i,:).*conj(dph_smooth(1:uw.n_ifg)));
    
    %dph_smooth_uw(1)=angle(dph_smooth(1));
    %for i2=2:uw.n_ifg
    %    dph_smooth_uw(i2)=dph_smooth_uw(i2-1)+angle(dph_smooth(i2).*conj(dph_smooth(i2-1)));
    %end
    dph_smooth_uw=cumsum(angle(dph_smooth(1:uw.n_ifg).*[1,conj(dph_smooth(1:uw.n_ifg-1))]));

    dph_close_master=mean(dph_smooth_uw(close_master_ix));
    dph_smooth_uw=dph_smooth_uw-(dph_close_master-angle(exp(j*dph_close_master)));
    
    dph_space_uw(i,:)=dph_smooth_uw+angle(dph_space(i,:).*conj(exp(j*dph_smooth_uw)));
end

if master_day~=0
    fprintf('\n   ESTIMATES OF DIFF PHASE NOISE STD DEV\n')
    fprintf('   =====================================\n')
    for i=1:uw.n_ifg
        fprintf('   %s  %4.1f deg\n',datestr(master_day+uw.day(i)),std(dph_noise(:,i))*180/pi)
    end
end

save('uw_unwrap_time','dph_space','dph_space_uw','close_master_ix','dph_noise');
