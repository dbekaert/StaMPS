function [] = masterchoice(Tcrit,Bcrit,Dcrit)
%   Script to calculate expected rho based on baseline information
%
%   K.H. Spaans, September 2010
%
%   ======================================================================
%   09/2010 K.H.: initial coding
%   ======================================================================


% Used by master_select
data = dlmread('Date.txt');
data(:,2) = dlmread('Bperp.txt');
data(:,3) = dlmread('Dopp.txt');

date = data(:,1);
[date,sortix] = sort(date);
data(:,1:2) = data(sortix,1:2);

data(:,1) = datenum(num2str(data(:,1)),'yyyymmdd');

for i=1:length(data(:,1))
    data2(:,1) = data(:,1) - data(i,1);
    data2(:,2) = data(:,2) - data(i,2);
    data2(:,3) = data(:,3) - data(i,3);
    T = (1-abs(data2(:,1)/Tcrit));
    B = (1-abs(data2(:,2)/Bcrit));
    D = (1-abs(data2(:,3)/Dcrit));
    T(T < 0) = 0;
    B(B < 0) = 0;
    D(D < 0) = 0;
    rho(i) = sum(T.*B.*D);
end

[dummy,sortix] = sort(rho,'descend');

fid = fopen('rho.txt','w');
fprintf(fid,'\nEstimated rho values:\n\n')
fprintf(fid,'Date:\t\trho:\n')
for i = 1:length(sortix)
	fprintf(fid,'%1.8g\t%2.2f\n',date(sortix(i)),rho(sortix(i)))
end
