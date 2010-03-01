%COREG_INV invert coreg translations for coefficients
%
%    Andy Hooper, Aug 2005
%
%   ======================================================================
%   04/2008 AH: Updated for compatibility with matlab 2008a
%   11/2008 AH: Deltaline/pixel added for Doris 3.96 compatibility
%   ======================================================================

cpmname=dir('CPM_Data.*');

load coreg_parms
Lrange=coreg_parms(1)-1
Prange=coreg_parms(2)-1
N=coreg_parms(3)^2
%osf=coreg_parms(4)
osf=1
n_ifg=coreg_parms(5)

G=sparse(0,n_ifg*12);
d=zeros(0,1);

for i=1:length(cpmname);
    
    thisname=cpmname(i).name;
    CPM_Data=load(thisname);
    
    posL1=((CPM_Data(:,2)-1)*4/Lrange)-2;
    posP1=((CPM_Data(:,3)-1)*4/Prange)-2;
    offL=CPM_Data(:,4);
    offP=CPM_Data(:,5);
    posL2=((CPM_Data(:,2)+offL-1)*4/Lrange)-2;
    posP2=((CPM_Data(:,3)+offP-1)*4/Prange)-2;
    n_pos=size(posL2,1);
    corrf=CPM_Data(:,6);
    stdev=sqrt(3/2/N).*sqrt(1-corrf.^2)/pi./corrf*osf^(3/2);
    weighting=1./stdev;
    ifgs=sscanf(thisname,'CPM_Data.%d.%d');
    ifg1=ifgs(1);
    ifg2=ifgs(2);
    
    Gblock1=[ones(size(posL1)),posL1,posP1,posL1.^2,posL1.*posP1,posP1.^2];
    Gblock1=Gblock1.*repmat(weighting,1,6);
    Gblock2=[ones(size(posL2)),posL2,posP2,posL2.^2,posL2.*posP2,posP2.^2];
    Gblock2=Gblock2.*repmat(weighting,1,6);
    Gnew=sparse(n_pos*2,n_ifg*12);
    if ifg1~=0
       Gnew(1:n_pos,(ifg1-1)*12+1:(ifg1-1)*12+6)=Gblock1;
       Gnew(n_pos+1:n_pos*2,(ifg1-1)*12+7:ifg1*12)=Gblock1;
    end
    Gnew(1:n_pos,(ifg2-1)*12+1:(ifg2-1)*12+6)=-Gblock2;
    Gnew(n_pos+1:n_pos*2,(ifg2-1)*12+7:ifg2*12)=-Gblock2;
    
    G=[G;Gnew];
    d=[d;weighting.*offL;weighting.*offP];
    %d=[d;offL;offP];
end

% coeff_s gives mapping of slave to master w.r.t. slave position
coeff_s=G\d;

load slave_corners.txt
deltaL=zeros(4,size(slave_corners,1));
deltaP=zeros(4,size(slave_corners,1));
for i = 1:size(slave_corners,1)
    posL1=((slave_corners(i,1:2))*4/Lrange)-2;
    posL1=[posL1(1);posL1(1);posL1(2);posL1(2)];
    posP1=((slave_corners(i,3:4))*4/Prange)-2;
    posP1=[posP1(1);posP1(2);posP1(1);posP1(2)];
    G=[ones(4,1),posL1,posP1,posL1.^2,posL1.*posP1,posP1.^2];
    deltaL(:,i)=G*coeff_s((i-1)*12+1:(i-1)*12+6);
    deltaP(:,i)=G*coeff_s((i-1)*12+7:(i)*12);
end
corner_offsets=[deltaL(:)';deltaP(:)'];
corner_offsets=corner_offsets(:);
save('corner_offsets.txt','-ascii','corner_offsets')

% we want mapping of slave to master w.r.t. master position
l=[-2:0.05:2];
p=l;
[Ls,Ps]=meshgrid(l,p);
Ls=Ls(:);
Ps=Ps(:);
nsynth=length(Ls);

Gblock=[ones(nsynth,1),Ls,Ps,Ls.^2,Ls.*Ps,Ps.^2];
Gsynth=sparse(2*nsynth,12);
Gsynth(1:nsynth,1:6)=Gblock;
Gsynth(nsynth+1:end,7:12)=Gblock;
coeff_m=zeros(size(coeff_s));

for i=1:n_ifg
    ifg_coeff_s=coeff_s((i-1)*12+1:i*12);
    dsynth=Gsynth*ifg_coeff_s;
    Lm=Ls+dsynth(1:nsynth)*4/Lrange;
    Pm=Ps+dsynth(nsynth+1:end)*4/Prange;
    Gblock=[ones(nsynth,1),Lm,Pm,Lm.^2,Lm.*Pm,Pm.^2];
    Gsynth(1:nsynth,1:6)=Gblock;
    Gsynth(nsynth+1:end,7:12)=Gblock;
    coeff_m((i-1)*12+1:i*12)=Gsynth\-dsynth;
end

save('coreg_coeffs.txt','-ascii','coeff_m')
      
    
