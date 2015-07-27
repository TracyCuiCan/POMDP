% This function returns the value and policy of multi-dimension sequential
% hypothesis testing, the input should be given like following form
% n=3; dimension
% M=20; how many part you discretise
% T=3;  horizon
% P=[0.25 0.75; 0.6 0.4; 0.7 0.3]; 2 observation, 3 hypothesis
% C=[20 20 20 1];  C0 C1¡­¡­Cn cc

function [V,G]=Discretize1(n,M,T,P,C)
%% initialization
points=multidimension(n,M); % call function multidimension to return discretized grids
obs=size(P,2);
point_num=size(points,2);
trans=zeros(point_num,point_num);
PTS=zeros(n+1,n);
Lambda=zeros(n+1,1);
V=zeros(T,point_num);
VC=zeros(T,point_num);
G=zeros(T,point_num);

%% Calculating Transition Probility
p_temp=zeros(n,1);
dist=zeros(point_num,1);
for i=1:point_num
    for j=1:obs
        for k=1:n
            p_temp(k)=P(k,j)*points(k,i)/(P(:,j)'*points(:,i)); % using bayes rule to get p+
        end
        for ii=1:point_num
            dist(ii)=mydistance(p_temp,points(:,ii)); % call function mydistance to calculate the distance between p+ and all grids 
        end
        [value, index]=sort(dist); % find the n nearest ones
        for jj=1:n  % p_temp=lambda_1*point1+lambda_2*point2+¡­¡­lambda*pointn ==> p_temp=[point1 point2 ... pointn]*[lambda_1; lambda_2;...lambda_n] also lambda1+lambda2+...lambdan=1
            PTS(:,jj)=[points(:,index(jj));1];
        end
        Lambda=pinv(PTS)*[p_temp;1];
        for kk=1:n
            trans(i,index(kk))=trans(i,index(kk))+Lambda(kk)*P(:,j)'*points(:,i); % update transition probability
        end
%         a=index(1);  % illustration for three dimension
%         b=index(2);
%         c=index(3);
%         trans(i,a)=trans(i,a)+lambda_a*P(:,j)'*points(:,i); 
%         trans(i,b)=trans(i,b)+lambda_b*P(:,j)'*points(:,i);
%         trans(i,c)=trans(i,c)+lambda_c*P(:,j)'*points(:,i);
    end
end

%%  Calculating Value function 
% terminal step
Vn=zeros(1,n);
for i=1:point_num % each point
    for j=1:n
        Vn(j)=points(j,i)*C(j); % each hypothesis
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
