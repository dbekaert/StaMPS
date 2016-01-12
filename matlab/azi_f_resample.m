function []=azi_f_resample(infile,informat,width,masterprf,slaveprf,outfile)
% Resample the slave image to the PRF of the master image
%
%   Hannes Bathke, February 2014
%   


if nargin<6
    help azi_f_resample
    error('Not enough arguments');
end

if strcmpi(informat,'complex_short')
    format='int16';
else 
    format='float';
end


n_patch = 3000;
overlap = n_patch/5;

buffer = zeros(n_patch,width,'single');
%n_fft=2048;
%l1=(n_fft-n_patch)/2;

% determine decimation (K) and interpolationfactors (L) satisfying (K/L = Min/Mout)
% and built filter kernel
resfac = masterprf/slaveprf;
[L M] = rat(resfac);


%Np = 5; % order
%ffrac = fdesign.polysrc(L,M,'Fractional Delay','Np',Np);
fnyq = fdesign.rsrc(L,M);
hfrac = design(fnyq,'kaiserwin');

fid=fopen(infile);                      % input slave slc
fid2=fopen(outfile,'w');                % resampled slave slc
fid3 = fopen('linenum.res','w');        % file for new line number

i=0;
ini_overlap = 0;
res_overlap = ini_overlap;
n_pix = 0;

disp('Resampling the slave.slc in azimuth to the PRF of the master.slc');
fprintf(['slave PRF: ',num2str(slaveprf),'\n']);
fprintf(['master PRF: ',num2str(masterprf),'\n']);
fprintf(['old patch size ',num2str(n_patch),'\n']);

while ~feof(fid)
    patch=fread(fid,[width*2,n_patch],[format,'=>single'])';
    n_patch=size(patch,1);
    new_patch_size = ceil(resfac*(n_patch+ini_overlap));
    disp(['new patch size ',num2str(new_patch_size)]);
    %if i == 4
        %break
    %    keyboard;
    %end
    if n_patch>0
        % no zeropadding because this changes the resfac and resulting
        % sampling frequency!!!
        
        if i > 0
            ov_buffer = zeros(n_patch+ini_overlap,width);
            ov_buffer(1:overlap,:)=buffer(prev_patch_size-ini_overlap+1:prev_patch_size,:);
            clear buffer
            buffer = ov_buffer;
        end
        
        buffer(ini_overlap+1:ini_overlap+size(patch,1),:)=complex(patch(:,1:2:end),patch(:,2:2:end));      
        clear patch
        
        fb = zeros(new_patch_size,width);
        for k = 1:width
            if k/500 == round(k/500)
                disp(['column ',num2str(k),' of patch ',num2str(i+1),]);
            end
            fb(:,k) = filter(hfrac,buffer(:,k));
        end
       
        %clear fb bb
        ifb=fb;
        
        patcht=zeros(new_patch_size-res_overlap,width*2,'single');
        patcht(:,1:2:end)=real(ifb(1+res_overlap:new_patch_size,:));
        patcht(:,2:2:end)=imag(ifb(1+res_overlap:new_patch_size,:));
        
        if (0)
            if i == 2
                keyboard
            
            figure; imagesc(abs(ifb(:,1:1000)))
            caxis([0 1000])
            figure; imagesc(abs(ifb(1+res_overlap:new_patch_size,1:1000)))
            caxis([0 1000])
            end
        end
        
        clear ifb 
        fwrite(fid2,patcht','float');
        i=i+1;
        
        fprintf('%d patch(es) processed\n',i)
        
        %overlap to previous patch
        ini_overlap = overlap;
        res_overlap = ceil(resfac*overlap);
        prev_patch_size = size(buffer,1);
        
        n_pix = size(patcht,1) + n_pix;
        fprintf('%d line(s) processed\n',n_pix)
    end
end

%calculate new number of lines

outpix = n_pix;

fprintf('Output lines: %d \n',outpix);
fprintf(fid3,'Output_lines: %d \n',outpix);

fclose(fid);
fclose(fid2);
fclose(fid3);

