function []=splitaz(infile,width,fdc,prf,informat,outfile)
%SPLITAZ Split azimuth spectrum in 2 and interfere the 2 halves
%
%   Andy Hooper, April 2006
%   

if nargin<6
    help splitaz
    error('Not enough arguments');
end

if strcmpi(informat,'complex_short')
    readformat='short=>single';
else
    readformat='float=>single';
end
    
n_patch=3000;
n_fft=4096;
l1=(n_fft-n_patch)/2;
filt1=zeros(n_fft,width,'single');
fdc_frac=fdc/prf-floor(fdc/prf+0.5)
bin1=round(n_fft/2+1+fdc_frac*n_fft);
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
    patch=fread(fid,[width*2,n_patch],readformat)';
    if size(patch,1)>0
        %buffer(l1+1:l1+size(patch,1),:)=complex(patch(:,1:2:end),patch(:,2:2:end));
        bb=complex(patch(:,1:2:end),patch(:,2:2:end));
        bb=fft(bb,n_fft);
        bb=ifft(bb.*filt2).*conj(ifft(bb.*filt1));
        bb=bb(1:size(patch,1),:);
        %bb=bb./abs(bb).*abs(patch);
        bb(isnan(bb))=0;
        patch(:,1:2:end)=real(bb);
        patch(:,2:2:end)=imag(bb);
        fwrite(fid2,patch','float');
        i=i+1;
        fprintf('%d patch(es) processed\n',i)
        %buffer(1:l1,:)=buffer(size(patch,1)+1:l1+size(patch,1),:);
        %buffer(l1+1:end,:)=0;
    end
end

fclose(fid);
fclose(fid2);
 
