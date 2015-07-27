% This function have three inputs, p0 and p1 are the
% probability of observation given hypothesis H0 and H1, grids are the
% disred number of grids. The return is the transitional probablilty
% following the principle of 1st order discretization.

function trans=Discretize0(P,grids)
obs=size(P,2);
N=grids; %divide [0,1] into N parts, grids
trans=zeros(N,N);
% probability, zero-hold order
points=linspace(0,1,N+1);
for i=1:length(points)-1
    beta= (points(i)+points(i+1))/2;
    for j=1:obs
        p_temp=beta*P(1,j)/(beta*P(1,j)+(1-beta)*P(2,j));
        k=max(1,ceil(p_temp*N));
        trans(i,k)=trans(i,k)+P(1,j)*beta+P(2,j)*(1-beta);
    end
    trans(i,:)=trans(i,:)/sum(trans(i,:));
end
end