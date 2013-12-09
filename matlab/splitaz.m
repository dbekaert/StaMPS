function []=splitaz(infile,informat,width,fdc,prf,outfile,plotflag)
%SPLITAZ Split azimuth spectrum in 2 and interfere the 2 halves
%
%   Andy Hooper, April 2006
%
%   ====================================================================
%   09/2006 AH: create all workspace files directly
%   12/2013 JvO: fixed bug for ascending images: subband filters shifted
%                based on the bin position of the max and min of meanbb
%   ====================================================================

if nargin<6
    help splitaz
    error('Not enough arguments');
end

if nargin<7
    plotflag=0;
end

if strcmpi(informat,'complex_short')
    format='int16';
else 
    format='float';
end

%infile(isnan(infile))=0+0i;
plotflag=1;    
n_patch=1500;
n_fft=2048;
l1=(n_fft-n_patch)/2;
buffer=zeros(n_fft,width,'single');
%fdc_frac=fdc/prf-floor(fdc/prf+0.5)
%bin1=round(n_fft/2+1+fdc_frac*n_fft);
%bin1=ceil(fdc_frac*n_fft);
%filt1=buffer;
%filt1=zeros(n_fft,1,'single');
%if fdc<0
%    filt1(bin1:bin1+n_fft/2-1,:)=1;
%else
%    filt1(bin1:end,:)=1;
%    filt1(1:bin1-n_fft/2-1,:)=1;
%end
%filt2=fftshift(filt1);

fid=fopen(infile);
fid2=fopen(outfile,'w');
fid3=fopen('bandwidthgap.txt','w');
i=0;

while ~feof(fid)
    patch=fread(fid,[width*2,n_patch],[format,'=>single'])';
    n_patch=size(patch,1);
    if n_patch>0
        buffer(l1+1:l1+size(patch,1),:)=complex(patch(:,1:2:end),patch(:,2:2:end));
        clear patch
        bb=fft(buffer);
        meanbb=sum(abs(bb),2);
        meanbb=meanbb/max(meanbb);
	minbb=find(meanbb==min(meanbb))
	maxbb=find(meanbb==max(meanbb))
        meanbb(meanbb<max(meanbb/5))=inf;
	if i==0
            n_bw=ceil(sum(meanbb~=inf)/2); % half bandwidth as number of bins (set 1st patch only);
	    fprintf(fid3,'Bandwidth gap = %f Hz\n',(n_bw)*prf/n_fft)
        end
        bin1=find(diff(meanbb)==inf)+1;
        bin1=bin1(1)-n_bw;
        if bin1<1
            bin1=bin1+n_fft/2;
        end
        filt=gausswin(n_bw,1.5);
        filt1=zeros(n_fft,1,'single');
        if bin1+n_bw-1<=n_fft
            filt1(bin1:bin1+n_bw-1,:)=filt;
        else
            filt1(bin1:end,:)=filt(1:n_fft-bin1+1);
            filt1(1:n_bw-(n_fft-bin1+1),:)=filt(n_fft-bin1+2:end);
        end
        if minbb < maxbb
		filt1=circshift(filt1,n_fft/2); %for ascending, disable for descending
	end
	filt2=circshift(filt1,-n_bw);

        if plotflag~=0
            fig=figure('Visible','off');
	    plot(meanbb)
            hold on
            plot(filt1./meanbb,'r')
            plot(filt2./meanbb,'g')
	    hold off
	    print(fig,'-dpng',sprintf(['splitaz_',int2str(i+1),'_',infile,'.png']));
            pause(0.1)
        end
        
        bb=ifft(bb.*repmat(filt2./meanbb,1,width)).*conj(ifft(bb.*repmat(filt1./meanbb,1,width)));
        %bb=ifft(bb.*repmat(filt2,1,width)).*conj(ifft(bb.*repmat(filt1,1,width)));
        absbb=abs(bb);
        absbb(absbb==0)=1;
        bb=bb./absbb.*abs(buffer);
        clear absbb
        patch=zeros(n_patch,width*2,'single');
        patch(:,1:2:end)=real(bb(l1+1:l1+size(patch,1),:));
        patch(:,2:2:end)=imag(bb(l1+1:l1+size(patch,1),:));
        clear bb
        fwrite(fid2,patch','float');
        i=i+1;
        fprintf('%d patch(es) processed\n',i)
        buffer(1:l1,:)=buffer(size(patch,1)+1:l1+size(patch,1),:);
        buffer(l1+1:end,:)=0;
    end
end

fclose(fid);
fclose(fid2);
fclose(fid3);
pause; 
