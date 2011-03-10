function llh=local2llh(xy,origin)
%local2llh     llh=local2llh(xy,origin)
%
%Converts from local coorindates to longitude and latitude 
%given the [lon, lat] of an origin. 'origin' should be in 
%decimal degrees. Note that heights are ignored and that 
%xy is in km.  Output is [lon, lat, height] in decimal 
%degrees. This is an iterative solution for the inverse of 
%a polyconic projection.

%-------------------------------------------------------------
%   Record of revisions:
%
%   Date          Programmer            Description of Change
%   ====          ==========            =====================
%
%   Aug 23, 2001  Jessica Murray        Clarification to help.
%
%   Apr 4, 2001   Peter Cervelli        Added failsafe to avoid
%                                       infinite loop because of
%                                       covergence failure.
%   Sep 7, 2000   Peter Cervelli		Original Code
%
%-------------------------------------------------------------

%Set ellipsoid constants (WGS84)

   a=6378137.0;
   e=0.08209443794970;

%Convert to radians / meters

%   xy=xy*1000;
%   origin=origin*pi/180;
   xy=double(xy)*1000;
   origin=double(origin)*pi/180;

%Iterate to perform inverse projection

   M0=a*((1-e^2/4-3*e^4/64-5*e^6/256)*origin(2) - ...
        (3*e^2/8+3*e^4/32+45*e^6/1024)*sin(2*origin(2)) + ...
        (15*e^4/256 +45*e^6/1024)*sin(4*origin(2)) - ...
        (35*e^6/3072)*sin(6*origin(2)));

   z=xy(2,:)~=-M0;

   A=(M0+xy(2,z))/a;
   B=xy(1,z).^2./a^2+A.^2;

   llh(2,z)=A;

   delta=Inf;

   c=0;
   
   while max(abs(delta))>1e-8

      C=sqrt((1-e^2*sin(llh(2,z)).^2)).*tan(llh(2,z));

      M=a*((1-e^2/4-3*e^4/64-5*e^6/256)*llh(2,z) - ...
           (3*e^2/8+3*e^4/32+45*e^6/1024)*sin(2*llh(2,z)) + ...
           (15*e^4/256 +45*e^6/1024)*sin(4*llh(2,z)) - ...
           (35*e^6/3072)*sin(6*llh(2,z)));

      Mn=1-e^2/4-3*e^4/64-5*e^6/256 - ...
         -2*(3*e^2/8+3*e^4/32+45*e^6/1024)*cos(2*llh(2,z)) + ...
         4*(15*e^4/256 +45*e^6/1024)*cos(4*llh(2,z)) + ...
         -6*(35*e^6/3072)*cos(6*llh(2,z));

      Ma=M/a;
   
      delta=-(A.*(C.*Ma+1)-Ma-0.5*(Ma.^2+B).*C)./ ...
           (e^2*sin(2*llh(2,z)).*(Ma.^2+B-2*A.*Ma)./(4*C)+(A-Ma).*(C.*Mn-2./sin(2*llh(2,z)))-Mn);

      llh(2,z)=llh(2,z)+delta;

      c=c+1;
      if c>100
          error('Convergence failure.')
      end
   end

   llh(1,z)=(asin(xy(1,z).*C/a))./sin(llh(2,z))+origin(1);

%Handle special case of latitude = 0

   llh(1,~z)=xy(1,~z)/a+origin(1);
   llh(2,~z)=0;

%Convert back to decimal degrees

   llh=llh*180/pi;
