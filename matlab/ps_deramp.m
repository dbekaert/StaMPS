function [ph_all,ph_ramp] = ps_deramp(ps,ph_all)
% [ph_all] = ps_deramp(ps,ph_all)
% Deramps the inputted data and gives that as output. Needs ps struct information!
%
% By David Bekaert
% modifications
% 09/2013   DB  Exclude nan-values in the plane estimation
% 09/2013   DB  Include the option when all points are NaN
% 11/2013   DB  Fix such deramped SM from SB inversion works as well
% 12/2013   DB  Add the ramp as output as well
% 05/2015   DB  Remove warning outputed to user
% 04/2016   DB  Make sure to use double to avoid machine precision warnings

fprintf(['Deramping computed on the fly. \n'])

% SM from SB inversion deramping
if ps.n_ifg ~= size(ph_all,2)
    ps.n_ifg = size(ph_all,2);
end

% detrenting of the data
A = [ps.xy(:,2:3)/1000 ones([ps.n_ps 1])];
ph_ramp = NaN(size(ph_all));
for k=1:ps.n_ifg
    ix = isnan(ph_all(:,k));
    if ps.n_ps-sum(ix)>5
        coeff = lscov(double(A(~ix,:)),double(ph_all(~ix,k)));
        ph_ramp(:,k)=A*coeff;
        ph_all(:,k)=ph_all(:,k)-ph_ramp(:,k);
    else
       fprintf(['Ifg ' num2str(k) ' is not deramped \n']) 
    end
end       