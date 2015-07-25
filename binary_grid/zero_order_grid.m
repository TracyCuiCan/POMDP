% This function finds the threshold for infinite horizon setup using zeroth-order discretization, 
% the discount parameters are beta = 1 and epsilon = 0.001;
tic

P=[0.25 0.75; 0.6 0.4];
obs=size(P,2);
cc=1;%cost to continue
L0=20; % the cost of testing H0 when H1 is true
L1=20; % the cost of testing H1 when H0 is true
N=1000; %divide [0,1] into N parts, grids
%M=1000; %discretize to approximate
trans=zeros(N,N);
% V=zeros(T,N);
% VC=zeros(T,N);
% G=zeros(T,N);

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

for i=1:N
    p(i)=0.5/N+(i-1)/N;
    V1=(1-p(i))*L1;
    V0=p(i)*L0;
    [V_zero(1,i),G(1,i)]=min([V0,V1]);
end

for n=2:1000
    for i=1:N
        Q=trans(i,:)*V_zero(n - 1,:)';
        VC(n,i)=cc+Q;
        [V_zero(n,i),G(n,i)]=min([p(i)*L0,(1-p(i))*L1,VC(n,i)]);
    end
    if max(abs(V_zero(n,:)-V_zero(n-1,:)))< 0.001 % check if the differnece is less than epsilon
        break
    end
    if n == 1000
        disp('Maximum iteration reached!')
    end
end
T = size(V_zero,1);
value = mean(V_zero(T,:));
lower = zeros(1, T);
upper = zeros(1, T);

lower(1)=(find(G(1,:)==2,1)-1)/N;
upper(1)=(find(G(1,:)==2,1)-1)/N;
for t=2:T
    lower(t)=(find(G(t,:)==3,1)-1)/N;
    upper(t)=(find(G(t,:)==2,1)-2)/N;
end
lower
upper

% figure
% for i=1
%  plot(p,VC(i,:),'b')
% % hold on
% %plot(p,V(i,:),'b','Linewidth',2)
% hold on
% xlabel('Grids');
% ylabel('Value Function');
% end
toc

            



