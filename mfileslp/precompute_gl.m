function P=precompute_gl(N,thsp_in,s_in,add,radius,Nslices)
Nspan=3;
beta=pi/Nspan;
dtheta_in=(thsp_in(2)-thsp_in(1));
ds_in=(s_in(2)-s_in(1));

%log-polar space
[Nrho,Ntheta,dtheta,drho,aR,am,g]=getparameters(beta,dtheta_in,ds_in,Nspan,N);
thsp=(-Ntheta/2:Ntheta/2-1)*dtheta;
rhosp=(-Nrho:-1)*drho;

thsp_in=thsp_in-beta/2;
for k=1:Nspan;
    pids{k}=find((thsp_in>=(k-1)*beta-beta/2).*(thsp_in<=(k-1)*beta+beta/2+1e-12));%get angles inside the current span
end
mids=Ntheta/2+1+(-ceil(Ntheta/4)-add:ceil(Ntheta/4)+add);
    

B3th=splineB3(thsp,radius);
B3rho=splineB3(rhosp,radius);
B3th=(fft(ifftshift(B3th)));
B3rho=(fft(ifftshift(B3rho)));
%save to structure
P=[];P.Nspan=Nspan;P.dtheta=dtheta;P.drho=drho;P.aR=aR;P.N=N;P.beta=beta;P.add=add;
P.thsp_in=thsp_in;P.s_in=s_in;P.Ns_in=length(s_in);P.Ntheta_in=length(thsp_in);P.pids=pids;
%theta
P.Ntheta=Ntheta;
P.Nthsp=numel(mids);
P.thsp0=thsp(1);
P.thspl=thsp(end)-thsp(1);
%rho
P.Nrho=Nrho;
P.rho=repmat(rhosp,numel(thsp),1);
P.rho0=rhosp(1);P.rhol=rhosp(end)-rhosp(1);
P.B3th=B3th;
P.B3rho=B3rho;
P.am=am;
P.g=g;
% P.Nslices=Nslices;