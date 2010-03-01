function calc_dem_offset()
%CALC_DEM_OFFSET calculate offset of dem
%
%   Andy Hooper, Dec 2007, modified after a script by Judicael Decriem
%                                             
%   ======================================================================
%   06/2009 AH: Large variables cleared from memory
%   10/2009 AH: More large variables cleared from memory
%   ======================================================================


im=combine_amp_dem;
im1=im(:,:,3);
im2=im(:,:,1);
clear im

[s1,s2] = size(im1);
im1 = single(im1);
im2 = single(im2);

a = fft2(im1-mean(mean(im1)));
clear im1  
b = ifft2(im2-mean(mean(im2)));
clear im2 
c = a.*b;
clear a b
c = ifft2(c);

[C,I1] = max(abs(c));
[C,I2] = max(max(abs(c)));

clear c

if(I2 > s2/2)
   dl = s2-I2+1;
else
   dl = 1-I2;
end   
if(I1(I2) > s1/2) 
    ds = s1-I1(I2)+1;
else
    ds = 1-I1(I2);
end    

fprintf('DEM down  = %d\n',ds)
fprintf('DEM right = %d\n',dl)
figure
plot_amp_dem(ds,dl);

