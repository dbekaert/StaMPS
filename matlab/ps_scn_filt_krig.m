function []=ps_scn_filt()
%PS_SCN_FILT estimate spatially correlated noise in unwrapped phase
%
%   Andy Hooper, June 2006
%
%   =======================================================================
%   11/2006 AH: Error corrected that was leaving master in temporal smoothing
%   04/2007 AH: Added 64-bit machine compatibility
%   05/2007 AH: Spatially correlated look angle error added
%   02/2010 AH: Replace unwrap_ifg_index with drop_ifg_index
%   02/2010 AH: Bug fixed that could cause a slave scn to be set to zero
%   11/2010 AH: If ramps estimated in step 7, subtract before scn estimation
%   02/2011 DB: Decreased time required for spatial filtering by factor 10
%   10/2011 MCC: Quadratic and seasonal defo is estimated per arc (edge) and then removed
%                before scn estimation (design matrix, A_temp, is hardcoded).
%                If eruption or earthquake a break point model is advicesable, but A_temp should be manually changed.
%
%   10/2011 MCC: Code to estimate scn_time_win and apply krigging for spatial estimation of APS
%   10/2011 MCC: Detects phase unwrapping errors per arc based on model and filtered defo.
%                Results are saved in unwrapping_errors_edges.mat
%   11/2011 MCC: Bug fixed the (probable) detected unwrapping errors are not removed from the original time series
%                Only saved for information
%   =======================================================================
logit;
fprintf('Estimating other spatially-correlated noise...\n')

pix_size=getparm('unwrap_grid_size',1);
time_win=getparm('scn_time_win',1);
deramp_ifg=getparm('scn_deramp_ifg',1);
scn_wavelength=getparm('scn_wavelength',1);
drop_ifg_index=getparm('drop_ifg_index',1);
small_baseline_flag=getparm('small_baseline_flag',1);
%krig_atmo=getparm('krig_atmo',1);
krig_atmo=true;

load psver
psname=['ps',num2str(psver)];
phuwname=['phuw',num2str(psver)];
sclaname=['scla',num2str(psver)];
apsname=['aps',num2str(psver)];
scnname=['scn',num2str(psver)]; % spatially-correlated noise
pmname=['pm',num2str(psver)];

ps=load(psname);
uw=load(phuwname);
%aps=load(apsname);

if strcmpi(small_baseline_flag,'y')
    unwrap_ifg_index=[1:ps.n_image];
else
    unwrap_ifg_index=setdiff([1:ps.n_ifg],drop_ifg_index);
end

day=ps.day(unwrap_ifg_index);
master_ix=sum(ps.master_day>ps.day)+1;
n_ifg=length(unwrap_ifg_index);
n_ps=ps.n_ps;

ph_all=single(uw.ph_uw(:,unwrap_ifg_index));
if exist([sclaname,'.mat'],'file')
    scla=load(sclaname);
    ph_all=ph_all-single(scla.ph_scla(:,unwrap_ifg_index));
    ph_all=ph_all-repmat(single(scla.C_ps_uw),1,length(unwrap_ifg_index));
    if ~isempty(scla.ph_ramp)
        ph_all=ph_all-single(scla.ph_ramp(:,unwrap_ifg_index));
    end
end
ph_all(isnan(ph_all))=0;

fprintf('   Number of points per ifg: %d',n_ps);

nodename=['scnfilt.1.node'];
fid=fopen(nodename,'w');
fprintf(fid,'%d 2 0 0\n',n_ps);

for i=1:n_ps
    fprintf(fid,'%d %f %f\n',i,ps.xy(i,2),ps.xy(i,3));
end

fclose(fid);

system('triangle -e scnfilt.1.node > triangle_scn.log');

fid=fopen('scnfilt.2.edge','r');
header=str2num(fgetl(fid));
N=header(1);
edges_nz=zeros(N,4);
for i=1:N
    edges_nz(i,:)=str2num(fgetl(fid));
end
fclose(fid);

%%% deramp end ifgs (unlike aps, orbit errors not so random and end
%%% orbit errors can pass through the low-pass filter
if strcmpi(deramp_ifg,'all') || krig_atmo
    deramp_ifg=1:ps.n_ifg;
end
deramp_ifg=intersect(deramp_ifg,unwrap_ifg_index);
deramp_ix=zeros(size(deramp_ifg));
ph_ramp=zeros(n_ps,length(deramp_ifg));

if ~isempty(deramp_ifg)
    fprintf('   deramping selected ifgs...\n')
    G=double([ones(n_ps,1),ps.xy(:,2),ps.xy(:,3)]);
    %G=double([ones(n_ps,1),ps.xy(:,2)]); % range only
    
    for i=1:length(deramp_ifg)
        i3=find(unwrap_ifg_index==deramp_ifg(i));
        deramp_ix(i)=i3;
        d=(ph_all(:,i3));
        m=G\double(d(:));
        ph_this_ramp=G*m;
        ph_all(:,i3)=ph_all(:,i3)-ph_this_ramp; % subtract ramp
        ph_ramp(:,i)=ph_this_ramp;
    end
    save(scnname,'ph_ramp')
end


%%% smooth in time using gaussian moving window
isnanix=isnan(uw.ph_uw);
uw.ph_uw(isnanix)=0;
dph=ph_all(edges_nz(:,3),:)-ph_all(edges_nz(:,2),:);
x_edges=(ps.xy(edges_nz(:,3),2)+ps.xy(edges_nz(:,2),2))*0.5;
y_edges=(ps.xy(edges_nz(:,3),3)+ps.xy(edges_nz(:,2),3))*0.5;


%%%%%%%%%%%%%%%%
%MCC starts here
%with respect to an arbitrary date doesnt matter which
years=(day- datenum('01-01-2000'))/365.25 ;
%mean_years=mean(years);

%thie is used to remove deformation before filtering
%MCC
%         quad                seasonal                          Const
A_time=[ years.^2 years sin(2*pi*years) cos(2*pi*years)-1 ones(size(years))];

AA=A_time'*A_time;

%dph_orig=dph;

modeled_defo_edges=repmat(single(NaN),size(A_time));

%%%%%%%%%%%%%%%%%
%This is used to calculate temporal variagrams
a=repmat(single(NaN),size(dph,1),1);
c1=repmat(single(NaN),size(dph,1),1);
c2=repmat(single(NaN),size(dph,1),1);
%inital values
a0=0.25;%range years
b0=1;%not used
c10=1;%noise variance rad
c20=2;%sill rad
cv_model=2;%2 exp; 3 gaussian
Nlags=25;%number of lags used for estimating variogram from histogram
plot_vario='n';%to do plots set to y
decor_dist=a0*6;%max time distance to include for variagram estimation in years
max_ob_vario=15000;
%%%%%%%%%%%%%%%%%%

rand_int=unique(randi([1,size(dph,1)],1,max_ob_vario));
if length(rand_int)<4000
  rand_int=unique(randi([1,size(dph,1)],1,max_ob_vario));
end
%this loop calculate defo from temporal model and at the same time the varaiogram per edge (arc)
tic
parfor n=1:size(dph,1)
    y=dph(n,:)';
    
    %parameters of modeld defo
    xhat=AA\A_time'*y;
    %residuals
    ehat=y-A_time*xhat;
    
    modeled_defo_edges(n,:)=xhat;
    
    dph(n,:)=ehat;
    %keyboard
    %fit a varigram per edge
    if ~isempty(find(rand_int==n, 1))
        [a(n) b c1(n) c2(n) emp_vario_defo_ps ] = ps_fit_vario( years,zeros(size(years)),dph(n,:),cv_model,a0,b0,c10,c20,decor_dist,plot_vario,Nlags);
    end
end
disp('Temporal variogram estimation:')

toc
save('modeled_defo_edges.mat','modeled_defo_edges','A_time','a','c1','c2');

fprintf('   low-pass filtering pixel-pairs in time...\n')

mean_x=nanmean(x_edges);
mean_y=nanmean(y_edges);
mean_std=nanmean(nanstd(dph));


%for some edges the estimation is not reliable
%We select only those estimated variagrams that are resaonable
ind_in=a<nanmean(years)*3 & a>nanmean(years)/50 & c1>mean_std/50  & c1<4*mean_std & c2>mean_std/50 &c2<4*mean_std & ~isnan(a);
save('modeled_defo_edges.mat','modeled_defo_edges','A_time','a','c1','c2','years','mean_std');


%keyboard
%if the number of remaninig edges is high enough (e.g. 100) we continue here otherwise filtering using the input window

dph_lpt=zeros(size(dph));
n_edges=size(dph,1);
kriging_method=2;

unwrapping_errors=repmat(int8(0),size(dph));
if length(find(ind_in))>50
    
    %We assume that the variagram coefficients (range, noise var and sill ) are spatially correlated and the spatial correlation can be described
    %with a second deg polynomial (designe matrix is A_vario.
    A_vario=[(x_edges(ind_in)/mean_x).^2  (y_edges(ind_in)/mean_y).^2  x_edges(ind_in)/mean_x  y_edges(ind_in)/mean_y  ones(size(x_edges(ind_in),1),1)  ];
    AA_vario=A_vario'*A_vario;
    
    %this are the coefficients describing range, noise var and sill spatially
    %this means that given an edge location, we are able to estimate its temporal varigram( range, noise var and sill)
    ahat=double(AA_vario)\A_vario'*a(ind_in);%range
    c1hat=double(AA_vario)\A_vario'*c1(ind_in);%noise var
    c2hat=double(AA_vario)\A_vario'*c2(ind_in);%sill
    
    Nmax=ceil(mean(years)+2*std(years));
    
    corrected_dph=repmat(single(NaN),size(dph));
    tic
    parfor n=1:n_edges
        
        %design matrix defining the temporal variao from edge position
        A_vario=[(x_edges(n)/mean_x).^2  (y_edges(n)/mean_y).^2  x_edges(n)/mean_x  y_edges(n)/mean_y  1 ];
        
        current_ahat=A_vario*ahat;
        current_c2hat=A_vario*c2hat;
        current_c1hat=A_vario*c1hat;
        %cv_model_current_edge

        if current_ahat<=nanmean(years)/50
          current_ahat=nanmean(a(ind_in));
        end
        if current_c2hat<=mean_std/50
          current_c2hat=nanmean(c2(ind_in));
        end
        if current_c1hat<=mean_std/50
          current_c1hat=nanmean(c1(ind_in));
        end
        cv_model_all=[1 NaN current_c2hat;cv_model current_ahat current_c1hat];
        
        %max decorelation distance depends on range (cv_model_all(2,2)) in years
        dx_max=cv_model_all(2,2)*4;
        
        temp_krigged=NaN(size(years));
        %kriging the random function (defo) from the varigram given by cv_model_all
        
        %we first check for unwrapping errors
        for nn=1:length(years)
           
            ind_in_time=find(abs(years-years(nn))<dx_max & abs(years-years(nn))~=0);
            [temp_krigged(nn) var_dph_lt  min_dist_year mean_dist_year ]=...
                ps_kriging([years(ind_in_time) zeros(size(ind_in_time)) dph(n,ind_in_time)' ] , [years(nn) 0 ] , Nmax,dx_max,kriging_method,cv_model_all);
            
        end%nn=1:length(years)
                
        %2*pi differences between current and predicted valua are assumed to be caused by unwrapping error
        unwrapping_errors(n,:)=round( (temp_krigged'-dph(n,:))/2.01/pi);
       
        corrected_dph(n,:)= dph(n,:)+single(unwrapping_errors(n,:))*2*pi;
        
        temp_krigged=NaN(size(years));
        
        for nn=1:length(years)
           
            ind_in_time=find(abs(years-years(nn))<dx_max );
            [temp_krigged(nn) var_dph_lt  min_dist_year mean_dist_year ]=...
                ps_kriging([years(ind_in_time) zeros(size(ind_in_time)) dph(n,ind_in_time)' ] , [years(nn) 0 ] , Nmax,dx_max,kriging_method,cv_model_all);
            
            
        end%nn=1:length(years)
        
        dph_lpt(n,:)=temp_krigged;
        
    end%end n=1:n_edges
    
    save('unwrapping_errors_edges.mat','x_edges','unwrapping_errors','y_edges');
    save('temp_vario.mat','ind_in','a','c1','c2','years','dph_lpt','dph','corrected_dph');
    
    disp('Time filtering :')
    toc
    %dph=corrected_dph;
    %MCC time filter ends
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else%if not enought data
    
    if isempty(time_win)%if time_win is empty we try to estimate it anyway
        time_win=2*nanmean(a(ind_in))*365.25;%a is in years and the filter is in days
    end
        
    for i1=1:n_ifg
        
        time_diff_sq=(day(i1)-day)'.^2;
        weight_factor=exp(-time_diff_sq/2/time_win^2);
        weight_factor(master_ix)=0; % leave out master
        weight_factor=weight_factor/sum(weight_factor);
        dph_lpt(:,i1)=sum(dph.*repmat(weight_factor,n_edges,1),2);
    end
end%if length(find(ind_in))>50


dph_hpt=dph-dph_lpt;  % leaves master APS - slave APS - slave noise (+ residue master noise)

ph_hpt=zeros(n_ps-1,n_ifg);
ref_ix=1;

A=sparse([[1:n_edges]';[1:n_edges]'],[edges_nz(:,2);edges_nz(:,3)],[-ones(n_edges,1);ones(n_edges,1)]);
A=double(A(:,[1:ref_ix-1,ref_ix+1:n_ps]));
%keyboard
fprintf('   solving for high-frequency (in time) pixel phase...\n')

for i=1:n_ifg
    ph_hpt(:,i)=A\double(dph_hpt(:,i));
end

ph_hpt=[ph_hpt(1:ref_ix-1,:);zeros(1,n_ifg);ph_hpt(ref_ix:end,:)]; % add back ref point

ph_scn=nan(n_ps,n_ifg);


ph_hpt=single(ph_hpt);

%%%%%%%%%%%%%%%%
%MCC starts atmo

if krig_atmo
    %coh_ps=getfield(load(pmname),'coh_ps');
    tic
    fprintf('   Kriging APS ...\n')
    
    a0=15000;%init range meters
    b0=1;%not used
    c10=1;%noise variance rad
    c20=1;%sill rad
    cv_model=2;%2 exp; 3 gaussian
    Nlags=30;%number of lags used for estimating variogram from histogram
    %plot_vario='n';%to do plots set to y
    decor_dist=a0*2;%max tim
    Nmax=50;
   
    %we cannot use all PS to estimate the variogram due to computer load
    %we select 10000 randomly distributed PS. Arbitrary decision but it should not matter
    % dist_to_mid=sqrt(( ps.xy(:,2)-mean_x_ps).^2 + ( ps.xy(:,3)-mean_y_ps).^2);
    % [sorted_dist ind_sort]=sort(dist_to_mid,'ascend');
    ind_vario= unique(randi(n_ps,min([max_ob_vario,n_ps]),1,'uint32'));
    if length(ind_vario)<4000 
				    ind_vario= unique(randi(n_ps,min([max_ob_vario,n_ps]),1,'uint32'));

    end
    %  ind_nearby_all=repmat(uint32(0),n_ps,Nmax+1);
%    orig_ind_vario=ind_vario;
    for n=1:n_ifg
 %       indnan=find(isnan(ph_hpt(ind_vario,n)));
        mean_x=nanmean(ps.xy(ind_vario,2));
        mean_y=nanmean(ps.xy(ind_vario,3));
        A_vario_ifg=[(ps.xy(ind_vario,2)/mean_x).^2  ps.xy(ind_vario,2)/mean_x  (ps.xy(ind_vario,3)/mean_y).^2 ps.xy(ind_vario,3)/mean_y ones(size(ps.xy(ind_vario,2)))];
        
        AA=double(A_vario_ifg'*A_vario_ifg);
        xhat=AA\A_vario_ifg'*ph_hpt(ind_vario,n);
        A_all=[(ps.xy(:,2)/mean_x).^2  ps.xy(:,2)/mean_x  (ps.xy(:,3)/mean_y).^2 ps.xy(:,3)/mean_y ones(size(ps.xy(:,2)))];
        %removes 2nd degree polynomial deterministic signal
        deramped_ph_hpt=ph_hpt(:,n)-A_all*xhat;
        
        [a b c1 c2 emp_vario_defo_ps ] = ps_fit_vario( ps.xy(ind_vario,2),  ps.xy(ind_vario,3),deramped_ph_hpt(ind_vario),cv_model,a0,b0,c10,c20,decor_dist,'n',Nlags);
        
        cv_model_all=[1 NaN c2;cv_model a c1];
        
        
        %max decorelation distance depends on range (cv_model_all(2,2)) in meters
        dx_max=cv_model_all(2,2)*4;
        tic
        aps=repmat(single(NaN),n_ps,1);
        %To improve performance I calculate the indeces of nearby PS only
        %once when n==1
        if n==1
            
            ind_nearby_all=repmat(uint32(0),n_ps,Nmax+1);
            for nn=1:n_ps
                dist=single(sqrt((ps.xy(:,2)-ps.xy(nn,2)).^2 + (ps.xy(:,3)-ps.xy(nn,3)).^2));
                [sorted_dist index]=sort(dist,'ascend');
                index2 = sorted_dist(1:Nmax)<dx_max;
                ind_nearby=  index(index2);
                ind_nearby_all(nn,1:length(ind_nearby))=ind_nearby;
                ind_nearby_all(nn,end)=length(ind_nearby);
                
                [aps(nn) var_dph_lt  min_dist_year mean_dist_year ]=...
                    ps_kriging([ ps.xy(ind_nearby,2) ps.xy(ind_nearby,3) deramped_ph_hpt(ind_nearby)], [ps.xy(nn,2) ps.xy(nn,3)] , Nmax,dx_max,kriging_method,cv_model_all);
            end%for length(years)
           
           
        else
            
            parfor nn=1:n_ps
                ind_nearby=ind_nearby_all(nn,1:ind_nearby_all(nn,end));
              
                [aps(nn) var_dph_lt  min_dist_year mean_dist_year ]=...
                    ps_kriging([ ps.xy(ind_nearby,2) ps.xy(ind_nearby,3) deramped_ph_hpt(ind_nearby)], [ps.xy(nn,2) ps.xy(nn,3)] , Nmax,dx_max,kriging_method,cv_model_all);
            end%parfor nn=1:length(years)
            
        end%end if n==1
        
        ph_scn(:,n)=aps+A_all*xhat;
        %figure;scatter(ps.xy(:,2), ps.xy(:,3),5,ph_scn(:,n),'filled');colorbar
        
        
    end%for n_ifg
    ph_scn=ph_scn+ph_ramp;
   
    %MCC APS ends
    toc
else% else krigged APS
    
    ph_hpt(:,deramp_ix)=ph_hpt(:,deramp_ix)+ph_ramp;
    
    ph_hpt=single(ph_hpt);
    
    sigma_sq_times_2=2*scn_wavelength.^2;
    patch_dist=scn_wavelength*4;
    patch_dist_sq=patch_dist*patch_dist;
    ix_range=ceil(n_ps/(max(ps.xy(:,3))-min(ps.xy(:,3)))*patch_dist*0.2);
    ix1=1;
    ix2=ix_range;
    ps.xy(:,1)=[1:n_ps]';
    
    fprintf('   low-pass filtering in space...\n')
    
    for i=1:n_ps
        
        x_min=ps.xy(i,2)-patch_dist;
        x_max=ps.xy(i,2)+patch_dist;
        y_min=ps.xy(i,3)-patch_dist;
        y_max=ps.xy(i,3)+patch_dist;
        
        ix1=ix1+ix_range;
        ix1(ix1>n_ps)=n_ps;
        while ix1>1 & ps.xy(ix1-1,3)>=y_min
            ix1=ix1-ix_range;
        end
        
        ix2=ix2-ix_range;
        ix2(ix2<1)=1;
        while ix2<n_ps & ps.xy(ix2+1,3)<=y_max
            ix2=ix2+ix_range;
        end
        
        ix1(ix1<1)=1;
        ix2(ix2>n_ps)=n_ps;
        
        xy_near=ps.xy(ix1:ix2,:);
        xy_near=xy_near(xy_near(:,2)>=x_min & xy_near(:,2)<=x_max & xy_near(:,3)>=y_min & xy_near(:,3)<=y_max,:);
        dist_sq=(xy_near(:,2)-ps.xy(i,2)).^2+(xy_near(:,3)-ps.xy(i,3)).^2;
        in_range_ix=dist_sq<patch_dist_sq; % exclude those out of range
        xy_near=xy_near(in_range_ix);
        dist_sq=dist_sq(in_range_ix);
        weight_factor=exp(-dist_sq/sigma_sq_times_2);
        weight_factor=weight_factor/sum(weight_factor); % normalize
        
        %ph_scn(i,:)=sum(ph_hpt(xy_near(:,1),:).*repmat(weight_factor,1,n_ifg),1); 		% slow method
        ph_scn(i,:)=weight_factor'*ph_hpt(xy_near(:,1),:);  				% speed increased by a factor 10
        
        if i/1000==floor(i/1000)
            disp([num2str(i),' PS processed'])
        end%if i/1000==floor(i/1000)
    end% end for i=1:n_ps
end % end if/else krig_atmo

ph_scn=ph_scn-repmat(ph_scn(1,:),n_ps,1); % re-ref to 1st PS
ph_scn_slave=zeros(size(uw.ph_uw));
ph_scn_slave(:,unwrap_ifg_index)=ph_scn;

ph_noise_res=zeros(size(uw.ph_uw));
ph_noise_res(:,unwrap_ifg_index)=ph_hpt-ph_scn_slave(:,unwrap_ifg_index);
std_ph_noise=nanstd(ph_noise_res(:,unwrap_ifg_index),0,2);
%figure;scatter(ps.xy(:,2), ps.xy(:,3),5,std_ph_noise,'filled');colorbar
%caxis([0 1])

ph_scn_slave(:,master_ix)=0;

save(scnname,'ph_scn_slave','ph_hpt','ph_ramp','ph_noise_res','std_ph_noise')
logit(1);
