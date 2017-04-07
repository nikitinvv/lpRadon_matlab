function Rsum_adj=Radjline(g,theta,s,N)
[xx1,xx2]=meshgrid(linspace(-1,1,N),linspace(-1,1,N));
Rsum_adj=zeros(N,N);
[T,S]=meshgrid(theta,s);
for j2=1:N
    for j1=1:N
            sp=xx1(j1,j2)*cos(theta)+xx2(j1,j2)*sin(theta);
            fp=qinterp2(T,S,g',theta,sp,2);           
            inds=find(isnan(fp)==0);
            Rsum_adj(j1,j2)=Rsum_adj(j1,j2)+sum(fp(inds)*(theta(2)-theta(1)));
    end
end
Rsum_adj=Rsum_adj*N;