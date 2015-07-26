%This function have three input, P is the conditionaly probability of
%observations given hypothesis H0 and H1; C is the cost matrix, C(1) is cost
%of deciding H1 while H0 is true, C(2) is cost of deciding H0 while H1 is
%true; C(3) is the continue cost
% The return is upper and lower bound of policy when the threshold remain
% unchanged

function Bound=find_bound(P,P_u,C)
global points
obs=size(P,2);
L1 = C(1);
L2 = C(2);
cc = C(3);
N = 10^3;
trans=zeros(N+1,N+1);
%points=linspace(0,1,N+1);

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
    V0= L1*P_u(i,1)*points(i) + L1*P_u(i,3)*(1-points(i)) + L2*P_u(i,2)*points(i);%u1 = 1
    V1= L1*P_u(i,2)*points(i) + L2*P_u(i,3)*(1-points(i)) + L1*P_u(i,4)*(1-points(i)); % u1 = 0
    [V(1,i),G(1,i)]=min([V0,V1]);
end
    
for n=2:1000
    for i=1:N+1
        Q=trans(i,:)*V(n - 1,:)';
        V0= L1*P_u(i,1)*points(i) + L1*P_u(i,3)*(1-points(i)) + L2*P_u(i,2)*points(i);%u1 = 1
        V1= L1*P_u(i,2)*points(i) + L2*P_u(i,3)*(1-points(i)) + L1*P_u(i,4)*(1-points(i)); % u1 = 0
        [V(n,i),G(n,i)]=min([V0,V1,cc+Q]);        
    end
    if max(abs(V(n,:)-V(n-1,:)))< 10^-3 % check if the differnece is less than epsilon
        break
    end
    if n == 1000
        disp('Maximum iteration reached!')
    end
end


T = size(V,1);
lower = (find(G(T,:)==1,1,'last'))/N;
upper = (find(G(T,:)==2,1)-2)/N;
% if lower > upper
%     lower = (find(G(T,:)==2,1)-1)/N;
%     upper = (find(G(T,:)==2,1)-1)/N;
% end
Bound = [lower upper];

end

