% This function is to evaluate and plot the performance of a certain policy
% specified by bound for infinite-horizon setup. 
tic
bound = [0.2 0.8];
P=[0.25 0.75; 0.6 0.4];
obs=size(P,2);
cc=1;%cost to continue
L0=20; % the cost of testing H0 when H1 is true
L1=20; % the cost of testing H1 when H0 is true
N=1000; %divide [0,1] into N parts, grids
trans=zeros(N+1,N+1);
lower = bound(1);
upper = bound(2);

% first-hold order
points=linspace(0,1,N+1);
for i=1:length(points)
    beta=points(i);
        for j=1:obs
            p_temp=beta*P(1,j)/(beta*P(1,j)+(1-beta)*P(2,j));
            k=max(1,ceil(p_temp*N));
            lambda=k-N*p_temp;
            p=P(1,j)*beta+P(2,j)*(1-beta);
            trans(i,k)=trans(i,k)+p*lambda;
            trans(i,k+1)=trans(i,k+1)+p*(1-lambda);
        end
    trans(i,:)=trans(i,:)/sum(trans(i,:));
end

V(1,1:N+1) = 0;

for n = 2:1000
    for i = 1:N+1
        if i < lower*N
            V(n,i) = points(i)*L0;
        else if i > upper*N
                V(n,i) = (1-points(i))*L1;
            else
                V(n,i) = trans(i,:)*V(n - 1,:)' + cc;
            end
        end
    end
    if max(abs(V(n,:)-V(n-1,:)))< 0.001 % check if the differnece is less than epsilon
        break
    end
end
T = size(V,1);
figure(2)
plot(points,V(T,:))
