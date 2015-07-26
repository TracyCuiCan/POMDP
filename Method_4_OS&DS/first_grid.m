% This function returns the optimal threshold and performance of
% decomposable decentralized hypothesis testing using value iteration based
% on first-order discretization. 

function [Bound, Value] = first_grid()
global trans1 trans2 C % The transition matrices and C are made global becuase this usually runs after fmin2.m 
N=10^3; %divide [0,1] into N parts, grids
%obs=size(P,2);
L0 = C(1);
L1 = C(1);
cc = C(3);

%% Sensor 1
for i=1:N+1
    p(i)=(i-1)/N;
    V1=(1-p(i))*L1;
    V0=p(i)*L0;
    [V_1(1,i),G_1(1,i)]=min([V0,V1]);
end

for n=2:1000
    for i=1:N+1
        Q=trans1(i,:)*V_1(n - 1,:)';
        VC(n,i) = cc + Q;
        [V_1(n,i),G_1(n,i)]=min([p(i)*L0,(1-p(i))*L1,VC(n,i)]);
    end
    if max(abs(V_1(n,:)-V_1(n-1,:)))< 10^-3 % check if the differnece is less than epsilon
        break
    end
end
T_1 = size(V_1,1);

%% sensor 2
for i=1:N+1
    p(i)=(i-1)/N;
    V1=(1-p(i))*L1;
    V0=p(i)*L0;
    [V_2(1,i),G_2(1,i)]=min([V0,V1]);
end

for n=2:1000
    for i=1:N+1
        Q=trans2(i,:)*V_2(n - 1,:)';
        VC(n,i) = cc + Q;
        [V_2(n,i),G_2(n,i)]=min([p(i)*L0,(1-p(i))*L1,VC(n,i)]);
    end
    if max(abs(V_2(n,:)-V_2(n-1,:)))< 10^-3 % check if the differnece is less than epsilon
        break
    end
    if n == 1000
        disp('Maximum iteration reached!')
    end
end
T_2 = size(V_2,1);

%% bound value
Bound = [(find(G_1(T_1,:)==3,1)-1)/N, (find(G_1(T_1,:)==2,1)-2)/N, (find(G_2(T_2,:)==3,1)-1)/N, (find(G_2(T_2,:)==2,1)-2)/N];
Value = mean(V_1(T_1,:),2) + mean(V_2(T_2,:),2);
end