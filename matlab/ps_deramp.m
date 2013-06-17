function [ph_all] = ps_deramp(ps,ph_all)

fprintf(['NOTE: Deramping computed on the fly. \n This will be changed latter on when tropospheric correction is incorporated in step 6 \n'])
% detrenting of the data
A = [ps.xy(:,2:3)/1000 ones([ps.n_ps 1])];
for k=1:ps.n_ifg
   coeff = lscov(A,ph_all(:,k));
   ph_all(:,k)=ph_all(:,k)-A*coeff;
end       