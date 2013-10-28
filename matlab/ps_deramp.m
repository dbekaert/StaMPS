function [ph_all] = ps_deramp(ps,ph_all)
% By David Bekaert
% modifications
% 09/2013   Exclude nan-values in the plane estimation
% 09/2013   Include the option when all points are NaN

fprintf(['NOTE: Deramping computed on the fly. \n This will be changed latter on when tropospheric correction is incorporated in step 6 \n'])
% detrenting of the data
A = [ps.xy(:,2:3)/1000 ones([ps.n_ps 1])];
for k=1:ps.n_ifg
    ix = isnan(ph_all(:,k));
    if ps.n_ps-sum(ix)>5
        coeff = lscov(A(~ix,:),ph_all(~ix,k));
        ph_all(:,k)=ph_all(:,k)-A*coeff;
    else
       fprintf(['Ifg ' num2str(k) ' is not deramped \n']) 
    end
end       