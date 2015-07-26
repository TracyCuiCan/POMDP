% This function finds the optimal thresholds for centralized hypothesis
% testing using policy iteration and the probability are computed using
% absoprtion probability. 
tic
P = [0.25 0.75; 0.6 0.4];
obs=size(P,2);
cc=1;%cost to continue
L0=20; % the cost of testing H0 when H1 is true
L1=20; % the cost of testing H1 when H0 is true
N=1000; %divide [0,1] into N parts, grids

global trans trans0 trans1
points = linspace(0,1,N+1);
trans=zeros(N+1,N+1);
trans0=zeros(N+1,N+1);
trans1=zeros(N+1,N+1);
for i=1:length(points)
    beta=points(i);
        for j=1:size(P,2)
            p_temp=beta*P(1,j)/(beta*P(1,j)+(1-beta)*P(2,j)); % updating p
            k=max(1,ceil(p_temp*N)); % check which part the posterior lies in
            lambda=k-N*p_temp;
            % P(i,j)
            p=P(1,j)*beta+P(2,j)*(1-beta); % linear interpolation
            trans(i,k)=trans(i,k)+p*lambda;
            trans(i,k+1)=trans(i,k+1)+p*(1-lambda);
            % P(i,j|H = 0)
            p=P(1,j)*beta;
            trans0(i,k)=trans0(i,k)+p*lambda;
            trans0(i,k+1)=trans0(i,k+1)+p*(1-lambda);
            % P(i,j|H=1)
            p=P(2,j)*(1-beta); 
            trans1(i,k)=trans1(i,k)+p*lambda;
            trans1(i,k+1)=trans1(i,k+1)+p*(1-lambda);
        end
    trans(i,:)=trans(i,:)/sum(trans(i,:));
    trans0(i,:)=trans0(i,:)/sum(trans0(i,:));
    trans1(i,:)=trans1(i,:)/sum(trans1(i,:));
end
trans0(1,:) = 0; trans0(1,1) = 1;
trans1(N+1,:) = 0; trans1(N+1,N+1) = 1;
trans = sparse(trans);
trans0 = sparse(trans0);
trans1 = sparse(trans1);


bound(1,:) = [0.5 0.5];
V(1,:) = absorption_policy_iter(bound(1,:));
W = zeros(1,N+1);    
for n=2:1000
    for i=1:N+1
        Q=trans(i,:)*V(n - 1,:)';
        [W(n,i),G(n,i)]=min([points(i)*L0,(1-points(i))*L1,cc + Q]);
    end
    bound(n,:) = [(find(G(n,:)==3,1)-1)/N, (find(G(n,:)==2,1)-2)/N];
    V(n,:) = absorption_policy_iter(bound(n,:));
    if bound(n,:) == bound(n-1,:)
        break
    end
end
toc
