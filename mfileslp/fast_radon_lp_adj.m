function [f]=fast_radon_lp_adj(R,Pgl,Padj)

f=zeros(Pgl.N);
R=[fliplr(R(end-Pgl.add:end-1,:)); R; fliplr(R(2:Pgl.add+1,:))];
s_in=Pgl.s_in;
thsp_in=[Pgl.thsp_in(1)-fliplr((1:Pgl.add)*(Pgl.thsp_in(2)-Pgl.thsp_in(1))) Pgl.thsp_in Pgl.thsp_in(end)+(1:Pgl.add)*(Pgl.thsp_in(2)-Pgl.thsp_in(1))];
rhosp=linspace(Pgl.rho0,Pgl.rho0+Pgl.rhol,Pgl.Nrho);
thsp=linspace(Pgl.thsp0,Pgl.thsp0+Pgl.thspl,Pgl.Ntheta);
mids=(Pgl.Ntheta-Pgl.Nthsp+1)/2+1:(Pgl.Ntheta+Pgl.Nthsp+1)/2;

thsp_in=single(thsp_in);s_in=single(s_in);R=single(R);rhosp=single(rhosp);thsp=single(thsp);Padj.fZintel=single(Padj.fZintel);
for k=1:Pgl.Nspan
    Padj.lp2p1{k}=single(Padj.lp2p1{k});Padj.lp2p2{k}=single(Padj.lp2p2{k});
    Padj.lp2p1w{k}=single(Padj.lp2p1w{k});Padj.lp2p2w{k}=single(Padj.lp2p2w{k}); 
    Padj.C2lp1{k}=single(Padj.C2lp1{k});Padj.C2lp2{k}=single(Padj.C2lp2{k});     
end
for k=1:Pgl.Nspan
   Rl1=single(zeros(Pgl.Ntheta,Pgl.Nrho));fp=single(zeros(Pgl.N));
   Rl1(Padj.lpids)=interp2(s_in,thsp_in,R,Padj.lp2p2{k},Padj.lp2p1{k},'spline',0); %interp polar->logpolar no borders!!
   fl=real(ifft2(fft2(Rl1).*Padj.fZ));%convolution (via FFT)
   fp(Padj.cids)=interp2(rhosp,thsp(mids),fl(mids,:),Padj.C2lp2{k},Padj.C2lp1{k},'spline',0); %interp logpolar->cart %will be wrapping
   f=f+fp;
end
f=f*(Pgl.N+1)/Pgl.N;

end
