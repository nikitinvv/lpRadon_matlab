function [F,G] = us2eqB3d(x,f,N)
x=mod(x-floor(x)+0.5+1,1)-0.5;
xeq1=(-N(1)/2:N(1)/2-1)';
xeq2=(-N(2)/2:N(2)/2-1);

phi=@(x)((1/6*(2-abs(x)).^3).*(abs(x)<2).*(abs(x)>=1)+(2/3-1/2*abs(x).^2.*(2-abs(x))).*(abs(x)<1));%need [-N/2:N/2-1]
B3theta=phi(xeq1);B3rho=phi(xeq2);B3thrho=B3theta*B3rho;
fB3=fftshift(fft2(ifftshift(B3thrho)));

M1=2;M2=2;
G=zeros(N(1)+2*M1,N(2)+2*M2);
phia=zeros(2*M1,2*M2);
for k=1:size(x,1)
  ell1=(floor(N(1)*x(k,1))-M1+1:floor(N(1)*x(k,1))+M1)';
  ell2=(floor(N(2)*x(k,2))-M2+1:floor(N(2)*x(k,2))+M2)';  
  for i2=1:2*M2
      for i1=1:2*M1
        phia(i1,i2)=phi(ell1(i1)-x(k,1)*N(1))*phi(ell2(i2)-x(k,2)*N(2));
      end
  end
  
  G(M1+1+N(1)/2+ell1,M2+1+N(2)/2+ell2)=G(M1+1+N(1)/2+ell1,M2+1+N(2)/2+ell2)+f(k)*phia;    
end
%wrap
[idx,idy]=ndgrid(-M1:N(1)+M1-1,-M2:N(2)+M2-1);
idx0=mod(idx+N(1),N(1))+1;idy0=mod(idy+N(2),N(2))+1;
G=accumarray([idx0(:),idy0(:)],reshape(G,[1,numel(G)]));
F=fftshift(fft2(fftshift(G)))./fB3;