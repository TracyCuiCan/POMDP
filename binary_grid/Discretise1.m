% This function have two inputs, P is the probability matrix, grids are the
% desired number of grids you want to discretize. The output is the N+1 by N+1 transitional probablilty
% following the principle of first order discretization.

function trans=Discretise1(P,N)
obs=size(P,2); 
trans=zeros(N+1,N+1);
% probability, first-hold order
points=linspace(0,1,N+1);
for i=1:length(points)
    beta=points(i);
        for j=1:obs
            p_temp=beta*P(1,j)/(beta*P(1,j)+(1-beta)*P(2,j)); % updating p
            k=max(1,ceil(p_temp*N)); % check which part the posterior lies in
            lambda=k-N*p_temp; 
            p=P(1,j)*beta+P(2,j)*(1-beta); % linear interpolation
            trans(i,k)=trans(i,k)+p*lambda;
            trans(i,k+1)=trans(i,k+1)+p*(1-lambda);
        end
    trans(i,:)=trans(i,:)/sum(trans(i,:));
end