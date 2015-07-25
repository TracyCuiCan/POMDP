% This function is to generate alpha-vector for infinite horizon scenario
% One can adapt this to different case by changing parameters.
% Normally, this function will take a long time to run when epsilon is small.
clear;
clc;
tic;

P=[0.25 0.75; 0.6 0.4];
obs=size(P,1); % size of Y
cc=1;%cost to continue
L0=20; % the cost of declaring H0 when H1 is true
L1=20; % the cost of declaring H1 when H0 is true

N = 1000;
points=linspace(0,1,N+1);
% when t=1, we only have first two term
for t=1
    Avec=[0 L1; L0 0]; % this is actually value function
    %Number(t)=size(Avec,1);
    x0=intersec([0 L1],[1,0],[0 0],[1 L0]);
    y0=[0 L1]*[x0; 1-x0];
    for n = 1:N+1
        for i = 1:size(Avec,1)
            V(n, i) = (Avec(i,2) - Avec(i,1))*points(n) + Avec(i,1);
        end
        W(1,n) = min(V(n,:));
    end
    taolower(1) = L1/(L0 + L1); %Both thresholds are the crosspoint of two lines
    taoupper(1) = L1/(L0 + L1);
end

for t=2:100
    %find all alpha-vectors under different Ys
    for j=1:obs
        temp(:,1)=Avec(:,1)*P(1,j); % a*f0(y)
        temp(:,2)=Avec(:,2)*P(2,j); %b*f1(y)
        a(j,:)=reshape(temp',1,size(temp,1)*size(temp,2));
    end
    
    % find all possible summation of alpha-vectors under different
    % observation Ys
    
    A=(reshape(a(1,:),2,length(a(1,:))/2))';
    for j=2:obs
        B=(reshape(a(j,:),2,length(a(j,:))/2))';
        A=minksum(A,B);
    end
    C=unique(A,'rows');
    ind=~any(C,2);
    C(ind,:)=[];
    
    % find the cross points of those vectors
    d=node(C);
    d=unique(d);
    d=[0,sort(d),1];
    
    % find minimum corresponding alpha-vectors
    for i=1:length(d)-1
        dtemp=(d(i)+d(i+1))/2;
        Ctemp=C*[dtemp; 1-dtemp];
        [V,index]=min(Ctemp,[],1);
        avec(i,1)=C(index,:)*[0;1];
        avec(i,2)=C(index,:)*[1;0];
    end
    avec=unique(avec,'rows');
    ind=~any(avec,2);
    avec(ind,:)=[];
    avec=avec+cc;
    
    % find the upper and lower bound for each horizon
    n=size(avec,1);
    for i=1:n
        upper(i)=intersec([x0, y0],[1,0],[0,avec(i,1)],[1,avec(i,2)]);
        lower(i)=intersec([0, 0],[x0,y0],[0,avec(i,1)],[1,avec(i,2)]);
    end
    lower=lower(find(lower>0 & lower<1));% finding lower and upper bound based on the index 
    upper=upper(find(upper>0 & upper<1));
    taolower(t)=min(lower);
    taoupper(t)=max(upper);
    
    
    %transform Avec back into a row vector and clear storage
    Avec=[avec; 0 L1; L0 0];
    
    for n = 1:N + 1
        for i = 1:size(Avec,1)
            V(n, i) = (Avec(i,2) - Avec(i,1))*points(n) + Avec(i,1);
        end
        W(t,n) = min(V(n,:));
    end
    epsilon = 0.1; 
    if max(abs(W(t,:)-W(t-1,:)))< epsilon % check if the differnece is less than epsilon, L1 norm
        break
    end
    
    clear a d B C  avec temp
end

% T = size(W,1)
% value = mean(W(T,:))
% 
% taolower
% taoupper
% 




toc

