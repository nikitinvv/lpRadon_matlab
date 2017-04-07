%testing version on CPU
clear all;
N=2^7;
Ntheta=3/2*N;Ns=N;
th=linspace(0,pi,Ntheta+1);th=th(1:end-1);s=linspace(-1,1,Ns);%initial polar grid
%% Parameters
osfZ=8;
Pgl=precompute_gl(N,th,s,4,1);
Padj=precompute_adj(Pgl,osfZ);

%circle to cut
[x1,x2]=meshgrid(linspace(-1,1,Pgl.N),linspace(-1,1,Pgl.N));circ0=(sqrt(x1.^2+x2.^2)<1-4/Pgl.N)*1.0;

disp('init exact filtered data');
%filtered phantom
[f,ellipse]=phantom(N);filter_kind='shepp-logan';%ramp,shepp-logan,cosine,cosine2,hamming,hann
ff=apply_filter_2d_exact(f,filter_kind,ellipse);
ff=ff.*circ0;
%filtered Radon data
h=apply_filter_exact(Ntheta,Ns,filter_kind,ellipse);h=single(h);

%% rec
disp('reconstruction');
frec=fast_radon_lp_adj(h,Pgl,Padj);frec=frec.*circ0;

%% result
figure(1);imagesc([frec ff]);title('true and reconstruction');
figure(2);imagesc(frec-ff);title('difference');
norm(abs(frec-ff)/norm(abs(ff),'fro'),'fro')