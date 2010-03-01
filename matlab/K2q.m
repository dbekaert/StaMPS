function [q]=K2q(K,range_pixel,sensor)
%[q]=K2q(K)
%convert from K to topo error
if nargin<2
    range_pixel=2300;
end
if nargin<3
    sensor='E';
end

if sensor=='A'
    %ALOS
    slant_range=4.684;
    near_range=846014.31;
    h=691650;
    lambda= 0.236057;
else
    %ERS/ENVISAT
    slant_range=7.9048902811596;
    near_range=831265;
    h=794141;
    lambda =0.0565647;
end

rho = near_range+range_pixel*slant_range;

re=6345245;
%load parms

alpha=pi-acos((rho.^2+re^2-(re+h)^2)/2./rho/re);

q=K*lambda.*rho.*sin(alpha)/4/pi;
