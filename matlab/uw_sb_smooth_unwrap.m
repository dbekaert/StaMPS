function [dph_smooth_series,F,model,energy,count]=uw_sb_smooth_unwrap(bounds,OPTIONS,G,W,dph,x1,varargin)
%UW_SB_SMOOTH_UNWRAP Unwrap by fitting a smooth function
%
%   Andy Hooper, June 2007


     if isempty(OPTIONS)
          scale=4;
          runs=3;
          grid=4;
          ts=linspace(1.5,2.5,runs);
          matrix=0;
          newton=0;
          talk=1;
     else
          OPTIONS(8)=0;
          if OPTIONS(1)
               scale=OPTIONS(1);
          else
               scale=4;
          end

          if OPTIONS(2)
               runs=OPTIONS(2);
          else
               runs=3;
          end

          if OPTIONS(3)
               grid=OPTIONS(3);
          else
               grid=4;
          end

          if OPTIONS(4)
               ts=ones(runs,1)*OPTIONS(4);
          else
               ts=linspace(2,3,runs);
          end

          matrix=OPTIONS(5);
          newton=OPTIONS(6);
          talk=OPTIONS(7);
     end

%Check bounds to make sure they're ok

     if max(bounds(:,1)>bounds(:,2))
          error('All the values in the first column of bounds must be less than those in the second.');
     end

%Define constants

     p=size(bounds,1);

     count=zeros(runs,1);
     energy=Inf;
     vals=[2.^-(1:grid)];
     vals=[vals,0,-vals];
     delta=0.5*abs((bounds(:,1)-bounds(:,2)));

%Loop through runs

n=size(G,2);

for k=1:runs
c=0;

     bestmodel=rand(p,100).*((bounds(:,2)-bounds(:,1))*ones(1,100))+bounds(:,1)*ones(1,100);

      O=zeros(100,1);
      for e=1:100
           %O(e)=feval(FUN,bestmodel(:,e),varargin{:});
           m=bestmodel(:,e);
           step_hat=m(1)*x1+m(2)*n/2*sin(2*pi/n*x1-m(3))+m(4)*n/2*sin(4*pi/n*x1-m(5));
                     dph_hat=G*step_hat;
                     dph_r=(dph-dph_hat)/2/pi;
                     dph_r=W*abs(dph_r-round(dph_r))*2*pi;
                     O(e)=sum(dph_r)+sum(abs(diff(step_hat)))/5;
      end
     tc=log10(mean(O))-ts(k);
     [v,i]=min(O);
     bestmodel=bestmodel(:,i);

if talk
     fprintf('\n\nBeginning run #%02d. Critical temperature at %3.2f.\n',k,10^tc);
     fprintf('------------------------------------------------\n\n');
     fprintf('f-Calls\t\tTemperature\tMinimum f-Value\n')
     fprintf('------------------------------------------------\n');
end

%Create cooling schedule from critical temperature and scale

     x=scale*[1 2 4 6 10 6 4 2 1];
     t=sum(x);
     temp=logspace(tc+1,tc-1,9);
     T=zeros(t,1);
     C=1;
     for i=1:9
          for j=1:x(i)
               T(C)=temp(i);
               C=C+1;
          end
     end

%Begin annealing

for w=1:t

     temp=T(w);
     c=c+1;

     if talk
          if c/10==floor(c/10)
               fprintf('%7d\t\t%7.2f\t\t%7.2f\n',count(k),temp,min(energy(1:c-1,k)));
          end
     end
     %Visit each parameter

      for x=1:p

          if delta(x)

          %Evaluate objective function at each permissible value

               v=bestmodel(x)+vals*delta(x);
               v=v(find((v<=bounds(x,2))&(v>=bounds(x,1))));
               mm=bestmodel*ones(1,length(v));
               mm(x,:)=v;
               NM=size(mm,2);
               count(k)=count(k)+NM;
                O=zeros(NM,1);
                for e=1:NM
                     %O(e)=feval(FUN,modelmatrix(:,e),varargin{:});
                     %m=mm(:,e);
                     step_hat=mm(1,e)*x1+mm(2,e)*n/2*sin(2*pi/n*x1-mm(3,e))+mm(4,e)*n/2*sin(4*pi/n*x1-mm(5,e));
                     dph_hat=G*step_hat;
                     dph_r=(dph-dph_hat)/2/pi;
                     dph_r=abs(W*(dph_r-round(dph_r)))*2*pi;
                     O(e)=sum(dph_r)+sum(abs(diff(step_hat)))/5;
                end

          %Form exponential probability distribution

               [dist,nanflag]=MakePDF(temp,O);
               if nanflag~=0
                    for u=1:length(nanflag)
                         disp(['Warning: the cost function generated a NaN for the following model:']);
                         disp(mm(nanflag(u)))
                    end
               end

          %Sample from probability distribution

               s=find(cumsum(dist)>=rand);
               if isempty(s)
                  fprintf('oops\n')
                  keyboard
               end
               s=s(1);
               bestmodel(x,1)=mm(x,s);
               energy(c,k)=O(s);
               model(:,c)=bestmodel;

          end
     end
end


[F(k,1),i]=min(energy(:,k));
mhat(:,k)=model(:,i);

if newton
     if newton==1
          % mstar=fminsearch(FUN,mhat(:,k),[],[],varargin{:});
          disp('FMINS is no longer a matlab function')
 	  % mstar=fmins(FUN,mhat(:,k),[],[],varargin{:});  %05-11-24 EKD
          Ostar=feval(FUN,mstar,varargin{:});
          if Ostar<F(k,1)
               if mstar>bounds(:,1) & mstar<bounds(:,2)
                    mhat(:,k)=mstar;
                    F(k,1)=Ostar;
                    if talk
                         fprintf('\nSimplex method lowered cost and remained within constraints.\n\n')
                    end
               else
                    if talk
                         fprintf('\nSimplex method lowered cost but failed to remain within constraints.\n\n')
                    end
               end
          end
     elseif newton==2
          mstar=constr(FUN,mhat(:,k),[],bounds(:,1),bounds(:,2),[],varargin{:});
          Ostar=feval(FUN,mstar,varargin{:});
          if Ostar<F(k,1)
             mhat(:,k)=mstar;
             F(k,1)=Ostar;
             if talk
                fprintf('\nConstrained optimization lowered cost.\n\n')
             end
          end

     end
end

end

[F,i]=min(F);
m=mhat(:,i);
dph_smooth_series=m(1)*x1+m(2)*n/2*sin(2*pi/n*x1-m(3))+m(4)*n/2*sin(4*pi/n*x1-m(5));


function [pdf,nanflag]=MakePDF(temp,v)
%Forms exponential probability distribution given a temperature and vector of costs.
%Internal function for simulated annealing algorithm.

bad=find(isnan(v));
if isempty(bad)
     pdf=eprob(temp,v);
     nanflag=0;
else
     good=find(~isnan(v));
     w=v(good);
     pdf=eprob(temp,w);
     pdf(good)=pdf;
     pdf(bad)=0;
     nanflag=bad;
end

function [pdf]=eprob(temp,v)
%Scales cost vector and calculates exponential probability distribution.  The scaling
%is necessary to permit wide ranges in temperature.
%Internal function for simulated annealing algorithm.

toobig=708.3964185322641;
pdf=v/temp;
mpdf=max(pdf);

if mpdf>toobig
   scale=mpdf/toobig;
   pdf=exp(-pdf/scale);
   pdf=pdf/max(pdf);
   pdf=pdf.^scale;
else
   pdf=exp(-pdf);
   pdf=pdf/max(pdf);
end

pdf=pdf/sum(pdf);

