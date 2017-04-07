function [Nrho,Ntheta,dtheta,drho,aR,am,g]=getparameters(beta,dtheta,ds,Nspan,N)
aR=sin(beta/2)/(1+sin(beta/2));
am=(cos(beta/2)-sin(beta/2))/(1+sin(beta/2));

g=[-0.09355974992572867 %Pi/3
   -0.05662612888153893 %Pi/4
   -0.03810640849602421 %Pi/5
   ];
g=g(Nspan-2);

Ntheta=N; %Nspan=3
Nrho=2*N;
 
dtheta=(2*beta)/Ntheta;
drho=(g-log(am))/Nrho;

end
