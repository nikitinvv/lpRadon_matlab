function g=Rphantom(ellipse,thsp,Ns)

rhosp=linspace(-1,1,Ns);
% thsp=linspace(0,pi,Ntheta);
[theta,rho]=ndgrid(thsp-pi/2,rhosp);
g=zeros(size(theta));




%
%   This head phantom is the same as the Shepp-Logan except 
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%      
%         A    a     b    x0    y0    phi
%        ---------------------------------
% % ellipse = [ 
% %     1  1   1    0     0     0
% % %     1   .69   .92    0     0     0   
% % %         -.8  .6624 .8740   0  -.0184   0
% % % %         0.5  .2100 .5100  .0    0    -70
% % % %         -0.5  .07100 .07100  .3    0.1    0
% %         -.2  .1600 .4100 -.22    0     18
% %          .1  .2100 .2500   0    .35    0
% % %          .1  .0460 .0460   0    .1     0
% % %          .1  .0460 .0460   0   -.1     0
% % %          .1  .0460 .0230 -.08  -.605   0 
% % %          .1  .0230 .0230   0   -.606   0
% % %          .1  .0230 .0460  .06  -.605   0 
% % % 1  0.01 0.01 0    0     0
% %          ];
  
   
g1=zeros(size(theta));
for k = 1:size(ellipse,1)    
   asq = ellipse(k,2)^2;       % a^2
   bsq = ellipse(k,3)^2;       % b^2
   phi = ellipse(k,6)*pi/180;  % rotation angle in radians
   x0 = ellipse(k,4);          % x offset
   y0 = ellipse(k,5);          % y offset
   A = ellipse(k,1);           % Amplitude change for this ellipse
% %    keyboard
   rho0=sqrt(asq*(cos(theta-phi)).^2+bsq*(sin(theta-phi)).^2);
   idx=find(abs(rho-x0*cos(theta)-y0*sin(theta))<=rho0);
   
   g(idx)=g(idx)+2*A*sqrt(asq*bsq)*sqrt(rho0(idx).^2-(rho(idx)-x0*cos(theta(idx))-y0*sin(theta(idx))).^2)./rho0(idx).^2;
   
   idx1=find(abs(rho-x0*cos(theta)-y0*sin(theta))<rho0);
   g1(idx1)=g1(idx1)+2*A*sqrt(asq*bsq)*sqrt(rho0(idx1).^2-(rho(idx1)-x0*cos(theta(idx1))-y0*sin(theta(idx1))).^2)./rho0(idx1).^2;
% keyboard;
end
