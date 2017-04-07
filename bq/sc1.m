clear all;
addpath('mfileslp','mfilter','/home/vnikitin/Dropbox/My/usflt');
   ellipse = [%1   .05  .05      0    0   0
            1    .04   .04    -0.5  -0.5     0
             %1   .03   .03    -.3  -.3   0
          ];
N=2^8;Ntheta=3*N/2;
th=linspace(0,pi,Ntheta+1);th=th(1:end-1);s=linspace(-1,1,N);%initial polar grid
%filtered phantom
[f,ellipse]=phantom(N);filter_kind='hamming';%ramp,shepp-logan,cosine,cosine2,hamming,hann
ff=apply_filter_2d_exact(f,filter_kind,ellipse);
%filtered Radon data
h=apply_filter_exact(Ntheta,N,filter_kind,ellipse);
Pgl=precompute_gl(N,th,s);
Padj=precompute_adj(Pgl,8);%osfZ=8

[x1,x2]=ndgrid(linspace(-1,1,N),linspace(-1,1,N));
circ=x1.^2+x2.^2<1-4/N;
frec=fast_radon_lp_adj(h,Pgl,Padj).*circ;
frecl=Radjline(h,th,s,N).*circ;

 
figure;imagesc([abs(frec-frecl)]);colorbar;
