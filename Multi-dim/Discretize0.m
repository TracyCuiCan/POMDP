% This function returns the value and policy of multi-dimension sequential
% hypothesis testing, the input should be given like following form
% n=3; dimension
% M=20; how many part you discretise
% T=3;  horizon
% P=[0.25 0.75; 0.6 0.4; 0.7 0.3]; 2 observation, 3 dimension
% C=[20 20 20 1];  C0 C1¡­¡­Cn cc

function [V,G]=Discretize0(n,M,T,P,C)

%% Initialization
points=multidimension(n,M); % call function multidimension to return discretized grids
obs=size(P,2);
point_num=size(points,2);
trans=zeros(point_num,point_num);
V=zeros(T,point_num);
VC=zeros(T,point_num);
G=zeros(T,point_num);

%% Transition Probility
p_temp=zeros(n,1);
dist=zeros(point_num,1);
for i=1:point_num
    for j=1:obs
        for k=1:n
            p_temp(k)=P(k,j)*points(k,i)/(P(:,j)'*points(:,i)); % using bayes rule to get p+
        end
        for ii=1:point_num
            dist(ii)=mydistance(p_temp,points(:,ii)); % call function mydistance to calculate the distance between p+ and eveay grid 
        end
        [dis, k]=min(dist); % find the closest point
        trans(i,k)=trans(i,k)+P(:,j)'*points(:,i); % update transition probability
    end
end

%% Value function 
% terminal step
Vn=zeros(1,n);
for i=1:point_num
    for j=1:n
        Vn(j)=points(j,i)*C(j);
    end
    [V(T,i),G(T,i)]=min(Vn);
end

% t=1:T-1
cc=C(length(C));
for t=T-1:-1:1
    for i=1:point_num
        Q=trans(i,:)*V(t+1,:)';
        VC(t,i)=cc+Q;
        for j=1:n
        Vn(j)=points(j,i)*C(j);
        end
        [V(t,i),G(t,i)]=min([Vn VC(t,i)]);
    end
end

end
