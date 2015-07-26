% This function finds the threshold for infinite horizon setup using first-order discretization, 
% the discount parameters are beta = 1 and epsilon = 0.001;

tic

P=[0.25 0.75; 0.6 0.4];
obs=size(P,2);
cc=1/20;%cost to continue
L0=1; % the cost of testing H0 when H1 is true
L1=1; % the cost of testing H1 when H0 is true
N=1000; %divide [0,1] into N parts, grids
%M=1000; %discretize to approximate
trans=zeros(N+1,N+1);
% V=zeros(T,N+1);
% VC=zeros(T,N+1);
% G=zeros(T,N+1);

% first-hold order
points=linspace(0,1,N+1);
%points = [0,unique(rand(1,N-1)),1];
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


for i=1:N+1
    p(i)=(i-1)/N;
    V1=(1-p(i))*L1;
    V0=p(i)*L0;
    [V_first(1,i),G(1,i)]=min([V0,V1]);
end
    
for n=2:1000
    for i=1:N+1
        Q=trans(i,:)*V_first(n - 1,:)';
        VC(n,i) = cc + Q;
        [V_first(n,i),G(n,i)]=min([p(i)*L0,(1-p(i))*L1,VC(n,i)]);
    end
    if max(abs(V_first(n,:)-V_first(n-1,:)))< 0.001 % check if the differnece is less than epsilon
        break
    end
    if n == 1000
        disp('Maximum iteration reached!')
    end
end
T = size(V_first,1);
value = mean(V_first(T,:));
lower = zeros(1, T);
upper = zeros(1, T);

lower(1)=(find(G(1,:)==2,1)-1)/N;
upper(1)=(find(G(1,:)==2,1)-1)/N;
for t=2:T %closed interval [lower, upper]
    lower(t)=(find(G(t,:)==3,1)-1)/N;  % index - 1
    upper(t)=(find(G(t,:)==2,1)-2)/N;
end
lower
upper
% figure(1)
% plot(points, V(T,:),'r');
% hold on
% cc=hsv(12);
% figure
% for i=T-1:-1:1
% % plot(p,VC(i,:),'b')
% % hold on
% 
% plot(p,V(i,:),'color',cc(i,:));
% hold on
% xlabel('Grids');
% ylabel('Value Function');
% end
toc