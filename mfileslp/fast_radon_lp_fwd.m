function [Rf fl Rl]=fast_radon_lp_fwd(f,Pgl,Pfwd)
Rf=zeros(Pgl.Ntheta_in,Pgl.Ns_in);

%form grids
x1=linspace(-1,1,Pgl.N);x2=linspace(-1,1,Pgl.N);
rhosp=linspace(Pgl.rho0,Pgl.rho0+Pgl.rhol,Pgl.Nrho);
thsp=linspace(Pgl.thsp0,Pgl.thsp0+Pgl.thspl,Pgl.Ntheta);
mids=(Pgl.Ntheta-Pgl.Nthsp+1)/2+1:(Pgl.Ntheta+Pgl.Nthsp+1)/2;

for k=1:Pgl.Nspan,
    fl=zeros(Pgl.Ntheta,Pgl.Nrho);
    [t,r]=ndgrid(x2,x1); plot(r(:),t(:),'r.');hold on;plot(Pfwd.lp2C1{k}(:),Pfwd.lp2C2{k}(:),'.');hold off;
    fl(Pfwd.cids)=interp2(x1,x2,f,Pfwd.lp2C1{k},Pfwd.lp2C2{k},'spline');%interp cart ->log-polar(circle)
%     file=fopen('t1','rb');fl1=fread(file,[Pgl.Ntheta, zPgl.Nrho],'float');fclose(file);
%     keyboard
    Rl=real(ifft2(fft2(fl.*exp(Pgl.rho)).*Pfwd.fZ));%convolution (via FFT).*exp(P.rho))
%     keyboard
    [t,r]=ndgrid(thsp(mids),rhosp); plot(r(:),t(:),'r.');hold on;plot(Pfwd.p2lp2{k}(:),Pfwd.p2lp1{k}(:),'.');hold off;
    Rs1=interp2(rhosp,thsp(mids),Rl(mids,:),Pfwd.p2lp2{k},Pfwd.p2lp1{k},'spline',0); %interp log-polar->polar
    Rf(Pgl.pids{k},:)=Rs1;
end
end
