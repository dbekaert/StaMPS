function [ph_uw]=uw_unwrap_space(max_rm_fraction)
%UW_UNWRAP_SPACE unwrap in space using output from UW_UNWRAP_TIME
%   as initial solution. 
%
%   Andy Hooper, June 2006
fprintf('Unwrapping in space...\n')

bad_std_thresh=10; % anything above is probably just noise

if nargin < 1
    max_rm_fraction=0.001; % maximum fraction of outlier edges to remove at one time
end

tic
%warning('ignore this warning')

uw=load('uw_phase');
ut=load('uw_unwrap_time');

res_tol=1e-3;         % tolerance for dph residual
ref_ix=1;

A=sparse(uw.n_edge,uw.n_ps);
for i=1:uw.n_edge
    A(i,uw.edgs(i,2))=-1;
    A(i,uw.edgs(i,3))=1;
end

A=A(:,[1:ref_ix-1,ref_ix+1:uw.n_ps]);
A_unweighted=A;


loop_sw=1;
defo_std=std(ut.dph_noise,0,2);
defo_std(defo_std>1.3)=defo_std(defo_std>1.3)+exp(defo_std(defo_std>1.3)-1.3)*4; % crude adjustment to covert wrapped std to unwrapped std

good_ix=defo_std<bad_std_thresh;

ix=good_ix;
n_dph=sum(ix);

success_flag=0;
i_loop=0;

ph_uw_ifg=zeros(uw.n_ps-1,1);
ph_uw=zeros(uw.n_ps,uw.n_ifg);

for i=1:uw.n_ifg
    disp(['PROCESSING IFG: ',num2str(i),' of ',num2str(uw.n_ifg)])
    
    ifg_std=exp(2*abs(ut.dph_noise(:,i)));
    %merged_std(merged_std<defo_std)=defo_std(merged_std<defo_std); % use higher of individual ifg error or series std.
    %weighting=1./merged_std;

    %weighting=1./defo_std;
    weighting=1./ifg_std;
    

    A=spdiags(weighting,0, uw.n_edge, uw.n_edge)*A_unweighted;
    
    dph_use=ut.dph_space_uw(:,i).*weighting;
    weighting_use=weighting;

    dph_use=dph_use(ix,:);
    edges_use=uw.edgs(ix,:);
    weighting_use=weighting_use(ix); % 
    A=A(ix,:);
    
    success_flag=0;
    edges_use=uw.edgs;
    n_dph=length(dph_use);
    n_bad=ceil(length(dph_use)*max_rm_fraction); % set n_bad to maximum fraction of bad edges that can be removed

    while success_flag==0

        if sprank(A)>=size(A,2)
            ph_uw_ifg=A\dph_use; % least squares inversion 
            A_save=A;
            edges_save=edges_use;
            weighting_save=weighting_use;
            dph_save=dph_use;
            n_dph_save=n_dph;
            dph_hat=A_save*ph_uw_ifg;
            r=(dph_save-dph_hat);                             % residual 
            %disp(['IFG ',num2str(i),': n_dph = ',num2str(n_dph_save),', residual std = ',num2str(std(r))]);
            fprintf('IFG %d: n_dph = %d, residual std = %f\n',i,n_dph_save,std(r));
        else
            ph_uw_ifg=A_save\dph_save; % least squares inversion 
            n_bad=ceil(n_bad/10);
            %disp(['IFG ',num2str(i),': n_dph = ',num2str(n_dph),' backing up']);

        end

        if abs(r) < res_tol % converged on an unwrapped solution
            success_flag=1;
        else
            [r_sort,sort_ix]=sort(abs(r)); % sort residuals in ascending order
            good_ix=logical(ones(n_dph_save,1));   
            ps_edge_dropped=zeros(uw.n_ps,1);

            % drop edges with the worst residuals, but only drop max of 1 edge per ps 
            for i2=length(sort_ix):-1:length(sort_ix)-n_bad+1
                    bad_ix=(sort_ix(i2));
                    if ps_edge_dropped(edges_save(bad_ix,2:3))==0  % if edge not already dropped for either ps of current edge
                        good_ix(bad_ix)=0;                        %   drop edge
                        ps_edge_dropped(edges_save(bad_ix,2:3))=1; %   mark ps as having an edge dropped
                    end
            end

            edges_use=edges_save(good_ix,:);
            weighting_use=weighting_save(good_ix);
            dph_use=dph_save(good_ix,:);
            A=A_save(good_ix,:);
            n_dph=length(dph_use);

        end
    
    end
    ph_uw(:,i)=[ph_uw_ifg(1:ref_ix-1);0;ph_uw_ifg(ref_ix:end)]; 
    save('uw_phaseuw','ph_uw')

end

% add back ref point
%ph_uw=[ph_uw_ifg(1:ref_ix-1,:);zeros(1,uw.n_ifg);ph_uw_ifg(ref_ix:end,:)]; 
ref_ph=angle(uw.ph(ref_ix,:));
ph_uw=ph_uw+repmat(ref_ph,uw.n_ps,1);

save('uw_phaseuw','ph_uw')

