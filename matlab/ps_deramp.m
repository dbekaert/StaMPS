function [ph_all,ph_ramp] = ps_deramp(ps,ph_all,degree)
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
% 12/2016   DB  Make sure a double is used for precision issues
% 01/2017   DB  Also allow for a 2nd order plane correction
% 03/2017   DB  Allow for a deramp_degree mat file to be specified

fprintf(['Deramping computed on the fly. \n'])

if nargin<3
    if exist('deramp_degree.mat','file')==2
        try 
            load('deramp_degree.mat','degree')
            fprintf('Found deramp_degree.mat file will use that value to deramp\n')
        catch
            degree = 1;
        end
    else
        degree = 1;
    end
end

% SM from SB inversion deramping
if ps.n_ifg ~= size(ph_all,2)
    ps.n_ifg = size(ph_all,2);
end

% detrenting of the data
if degree == 1
    % z = ax + by+ c
    A = double([ps.xy(:,2:3)/1000 ones([ps.n_ps 1])]);    
    fprintf('**** z = ax + by+ c\n')
elseif degree == 1.5
    % z = ax + by+ cxy + d
    A = double([ps.xy(:,2:3)/1000 ps.xy(:,2)/1000.*ps.xy(:,3)/1000 ones([ps.n_ps 1])]);    
    fprintf('**** z = ax + by+ cxy + d\n')
elseif degree ==2
    % z = ax^2 + by^2 + cxy + d
    A = double([(ps.xy(:,2:3)/1000).^2 ps.xy(:,2)/1000.*ps.xy(:,3)/1000 ones([ps.n_ps 1])]);
    fprintf('**** z = ax^2 + by^2 + cxy + d \n')
elseif degree ==3
    % z = ax^3 + by^3 + cx^2y + dy^2x + ex^2 + by^2 + cxy + d
    A = double([(ps.xy(:,2:3)/1000).^3  (ps.xy(:,2)./1000).^2.*ps.xy(:,3)/1000 (ps.xy(:,3)./1000).^2.*ps.xy(:,2)/1000  ...
                (ps.xy(:,2:3)/1000).^2 ps.xy(:,2)/1000.*ps.xy(:,3)/1000 ones([ps.n_ps 1])]);
    fprintf('**** z = ax^3 + by^3 + cx^2y + dy^2x + ex^2 + fy^2 + gxy + h \n')
    
end

ph_ramp = NaN(size(ph_all));
for k=1:ps.n_ifg
    ix = isnan(ph_all(:,k));
    if ps.n_ps-sum(ix)>5
        coeff = lscov(A(~ix,:),ph_all(~ix,k));
        ph_ramp(:,k)=A*coeff;
        ph_all(:,k)=ph_all(:,k)-ph_ramp(:,k);
    else
       fprintf(['Ifg ' num2str(k) ' is not deramped \n']) 
    end
end       
