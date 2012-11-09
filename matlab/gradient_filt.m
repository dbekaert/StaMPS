function [ifreq,jfreq,ij,Hmag]=gradient_filt(ph,n_win)
%GRADIENT_FILT Determine 2-D gradient through FFT
%   [ph_out]=gradient_filt(ph,n_win)
%
%   Andy Hooper, Feb 2012
%   


[n_i,n_j]=size(ph);
n_inc=floor(n_win/4);
n_win_i=ceil(n_i/n_inc)-3;
n_win_j=ceil(n_j/n_inc)-3;

ph(isnan(ph))=0;
%B=gausswin(3)*gausswin(3)';
%B=B./sum(B(:));
ph_bit=zeros(n_win);

Hmag=nan(n_win_i,n_win_j);
ifreq=Hmag;
jfreq=Hmag;
ij=nan(n_win_i*n_win_j,2);
i=0;

for ix1=1:n_win_i
    i1=(ix1-1)*n_inc+1;
    i2=i1+n_win-1;
    if i2>n_i
        i_shift=i2-n_i;
        i2=n_i;
        i1=n_i-n_win+1;
    end
    for ix2=1:n_win_j
        i=i+1;
        j1=(ix2-1)*n_inc+1;
        j2=j1+n_win-1;
        if j2>n_j
           j_shift=j2-n_j;
           j2=n_j;
           j1=n_j-n_win+1;
        end
        ph_bit(1:n_win,1:n_win)=ph(i1:i2,j1:j2);
        if sum(ph_bit(:)~=0)<6 % cannot reliably estimate gradient from fewer points
            Hmag(ix1,ix2)=nan;
            ifreq(ix1,ix2)=nan;
            jfreq(ix1,ix2)=nan;
        else
            ph_fft=fft2(ph_bit);
            %H=abs(fftshift(ph_fft));
            %H=filter2(B,H); % smooth response  
            H=abs(ph_fft);
            [Hmag_this,I]=max(H(:));
            Hmag(ix1,ix2)=Hmag_this/mean(H(:));
            
            I1=mod(I,n_win);
            I2=ceil(I/n_win);
            
            I1=mod(I1+n_win/2,n_win);
            I1(I1==0)=n_win;
            I2=mod(I2+n_win/2,n_win);
            I2(I2==0)=n_win;
            
            ifreq(ix1,ix2)=(I1-n_win/2-1)*2*pi/n_win; % i is top to bottom, fft is bottom to top
            jfreq(ix1,ix2)=(I2-n_win/2-1)*-2*pi/n_win; 
        end
        ij(i,:)=[(i1+i2)/2,(j1+j2)/2];
    end
end
ifreq=ifreq';
jfreq=jfreq';
Hmag=Hmag';

