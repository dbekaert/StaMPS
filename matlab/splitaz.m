function []=splitaz(infile,width,fdc,prf,outfile)
%SPLITAZ Split azimuth spectrum in 2 and interfere the 2 halves
%
%   Andy Hooper, April 2006
%   


if nargin<5
    help splitaz
    error('Not enough arguments');
end

    
n_patch=3000;
n_fft=4096;
l1=(n_fft-n_patch)/2;
buffer=zeros(n_fft,width,'single');
fdc_frac=fdc/prf-floor(fdc/prf+0.5)
bin1=round(n_fft/2+1+fdc_frac*n_fft);
filt1=buffer;
if fdc<0
    filt1(bin1:bin1+n_fft/2-1,:)=1;
else
    filt1(bin1:end,:)=1;
    filt1(1:bin1-n_fft/2-1,:)=1;
end
filt2=fftshift(filt1);

fid=fopen(infile);
fid2=fopen(outfile,'w');
i=0;

while ~feof(fid)
    patch=fread(fid,[width*2,n_patch],'float')';
    if size(patch,1)>0
        buffer(l1+1:l1+size(patch,1),:)=complex(patch(:,1:2:end),patch(:,2:2:end));
        bb=fft(buffer);
        bb=ifft(bb.*filt2).*conj(ifft(bb.*filt1));
        bb=bb./abs(bb).*abs(buffer);
        patch(:,1:2:end)=real(bb(l1+1:l1+size(patch,1),:));
        patch(:,2:2:end)=imag(bb(l1+1:l1+size(patch,1),:));
        fwrite(fid2,patch','float');
        i=i+1;
        fprintf('%d patch(es) processed\n',i)
        buffer(1:l1,:)=buffer(size(patch,1)+1:l1+size(patch,1),:);
        buffer(l1+1:end,:)=0;
    end
end

fclose(fid);
fclose(fid2);
 
