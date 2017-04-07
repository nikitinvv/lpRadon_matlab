function Padj=precompute_adj(P,osfZ)
Nspan=P.Nspan;beta=P.beta;Ntheta=P.Ntheta;Nrho=P.Nrho;drho=P.drho;dtheta=P.dtheta;aR=P.aR;

thsp=(-Ntheta/2:Ntheta/2-1)*dtheta;
rhosp=(-Nrho:-1)*drho;
%filter
fZ=(fftshift(fzeta_loop_weights_adj(Ntheta,Nrho,Ntheta*dtheta,Nrho*drho,0,osfZ)))*2/(P.s_in(2)-P.s_in(1)); %last parameter - os for angles 
C2lp1=[];C2lp2=[];
[x1,x2]=meshgrid(linspace(-1,1,P.N),linspace(-1,1,P.N));
cids=find(x1.^2+x2.^2<=1);
for k=1:Nspan;
    z1=aR*( x1(cids)*cos((k-1)*beta+beta/2)+x2(cids)*sin((k-1)*beta+beta/2))+(1-aR);
    z2=aR*(-x1(cids)*sin((k-1)*beta+beta/2)+x2(cids)*cos((k-1)*beta+beta/2));
    C2lp1{k}=atan2(z2,z1);C2lp2{k}=log(sqrt(z1.^2+z2.^2));  
end

% keyboard;
[z1,z2]=ndgrid(thsp,exp(rhosp));

z2n=z2-(1-aR)*cos(z1);z2n=z2n/aR;
lpids=find((z1>=-beta/2).*(z1<beta/2).*(abs(z2n)<=1));%cut in theta and rho
% wids=r
lp2p1=[];lp2p2=[];
for k=1:Nspan;
lp2p1{k}=z1(lpids)+(k-1)*beta;
lp2p2{k}=z2n(lpids);
end

%right side
wids=find((log(z2)>+P.g));
z2n=exp(log(z2(wids))+log(P.am)-P.g)-(1-aR)*cos(z1(wids));z2n=z2n/aR;
lpidsw=find((z1(wids)>=-beta/2).*(z1(wids)<beta/2).*(abs(z2n)<=1));%cut in theta and rho

%left side
wids2=find(log(z2)<log(P.am)-P.g+drho);
z2n2=exp(log(z2(wids2))-log(P.am)+P.g)-(1-aR)*cos(z1(wids2));z2n2=z2n2/aR;
lpidsw2=find((z1(wids2)>=-beta/2).*(z1(wids2)<beta/2).*(abs(z2n2)<=1+P.add*drho));%cut in theta and rho

for k=1:Nspan;
lp2p1w{k}=z1([lpidsw;lpidsw2])+(k-1)*beta;
lp2p2w{k}=[z2n(lpidsw);z2n2(lpidsw2)];
end
for k=1:Nspan
    thsp_in0(k)=P.thsp_in(P.pids{k}(1));
    thsp_inl(k)=P.thsp_in(P.pids{k}(end))-P.thsp_in(P.pids{k}(1));
end

p2lp1=[];p2lp2=[];
for k=1:Nspan;
    thsp0=P.thsp_in(P.pids{k});
    [th0,s0]=ndgrid(thsp0-(k-1)*beta,P.s_in);
    p2lp1{k}=th0;
    p2lp2{k}=log(s0*aR+(1-aR)*cos(th0));
end


B3com=P.B3th'*P.B3rho;
%save to structure
P.adj=[];Padj.Nspan=Nspan;Padj.fZ=fZ;Padj.fZgpu=Padj.fZ(1:end/2+1,:)./B3com(1:end/2+1,:);
Padj.fZintel=Padj.fZ(1:end/2+1,:);Padj.fZintel=conj(Padj.fZintel);Padj.fZintel=[real(Padj.fZintel(:)');imag(Padj.fZintel(:)')];
Padj.lp2p1=lp2p1;Padj.lp2p2=lp2p2;
Padj.p2lp1=p2lp1;Padj.p2lp2=p2lp2;
Padj.lp2p1w=lp2p1w;Padj.lp2p2w=lp2p2w;
Padj.C2lp1=C2lp1;Padj.C2lp2=C2lp2;
Padj.cids=cids;
Padj.lpids=lpids;
Padj.thsp_in0=thsp_in0;
Padj.thsp_inl=thsp_inl;
Padj.s_in0=P.s_in(1);
Padj.s_inl=P.s_in(end)-P.s_in(1);
Padj.wids=[wids(lpidsw);wids2(lpidsw2)];
