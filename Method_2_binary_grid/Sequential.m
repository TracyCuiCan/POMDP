%This function have five input, P0 and P1 are conditionaly probability of
%observations given hypothesis H0 and H1; C is continuous cost, C0 is cost
%of deciding H1 while H0 is true, C1 is cost of deciding H0 while H1 is
%true; hoizon is the desired time step
% The return is upper and lower bound of policy

function [lower, upper]=Sequential(P,C,horizon)
T=horizon;
obs=size(P,2);
cc=C(3);%cost to continue
L0=C(1); % the cost of testing H0 when H1 is true
L1=C(2); % the cost of testing H1 when H0 is true
N=1000; %divide [0,1] into N parts, grids
M=1000; %discretize to approximate
trans=zeros(N,N);
V=zeros(T,N);
G=zeros(T,N);
lower=zeros(T,1);
upper=zeros(T,1);

% probability, zero-hold order
points=linspace(0,1,N+1);
for i=1:length(points)-1
    beta=linspace(points(i),points(i+1),M+1);
    p(i)=0.5/N+(i-1)/N;
    for beta=beta(2:length(beta)-1)
        for j=1:obs
            p_temp=beta*P(1,j)/(beta*P(1,j)+(1-beta)*P(2,j));
            k=max(1,ceil(p_temp*N));
            trans(i,k)=trans(i,k)+P(1,j)*p(i)+P(2,j)*(1-p(i));
        end
    end
    trans(i,:)=trans(i,:)/sum(trans(i,:));
end

for i=1:N
    p(i)=0.5/N+(i-1)/N;
    V1=(1-p(i))*L1;
    V0=p(i)*L0;
    [V(T,i),G(T,i)]=min([V0,V1]);
end
    
for n=T-1:-1:1
 for i=1:N
     Q=trans(i,:)*V(n+1,:)';
     [V(n,i),G(n,i)]=min([p(i)*L0,(1-p(i))*L1,cc+Q]);
 end
end

lower(1)=(find(G(T,:)==2,1)-1)/N;
upper(1)=(find(G(T,:)==2,1)-1)/N;



for t=1:T-1
    lower(T-t+1)=(find(G(t,:)==3,1)-1)/N;
    upper(T-t+1)=(find(G(t,:)==2,1)-1)/N;
end
