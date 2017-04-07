function Pfwd=precompute_fwd(P,osfZ)
Nspan=P.Nspan;beta=P.beta;Ntheta=P.Ntheta;Nrho=P.Nrho;drho=P.drho;dtheta=P.dtheta;aR=P.aR;

thsp=(-Ntheta/2:Ntheta/2-1)*dtheta;
rhosp=(-Nrho:-1)*drho;

%filter
fZ=fftshift(fzeta_loop_weights(Ntheta,Nrho,Ntheta*dtheta,Nrho*drho,0,osfZ))/aR;

lp2C1=[];lp2C2=[];
for k=1:Nspan;
  lp2C1{k}=((cos(thsp')*exp(rhosp)-(1-aR))*cos((k-1)*beta+beta/2)-sin(thsp')*exp(rhosp)*sin((k-1)*beta+beta/2))/aR;
  lp2C2{k}=((cos(thsp')*exp(rhosp)-(1-aR))*sin((k-1)*beta+beta/2)+sin(thsp')*exp(rhosp)*cos((k-1)*beta+beta/2))/aR;
  cids=find((lp2C1{k}.^2+lp2C2{k}.^2)<=1);
  lp2C1{k}=lp2C1{k}(cids);
  lp2C2{k}=lp2C2{k}(cids);
end

p2lp1=[];p2lp2=[];
for k=1:Nspan;
    thsp0=P.thsp_in(P.pids{k});
    [th0,s0]=ndgrid(thsp0-(k-1)*beta,P.s_in);
    p2lp1{k}=th0;
    p2lp2{k}=log(s0*aR+(1-aR)*cos(th0));
end

%save to structure
B3com=P.B3th'*P.B3rho;
Pfwd=[];Pfwd.Nspan=Nspan;Pfwd.fZ=fZ;Pfwd.fZgpu=Pfwd.fZ(1:end/2+1,:)./(B3com(1:end/2+1,:));%todo
Pfwd.lp2C1=lp2C1;Pfwd.lp2C2=lp2C2;
Pfwd.p2lp1=p2lp1;Pfwd.p2lp2=p2lp2;
Pfwd.cids=cids;
Pfwd.x0=-1;Pfwd.xl=2;
Pfwd.y0=-1;Pfwd.yl=2;