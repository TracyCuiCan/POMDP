% This gives optimal bounds for three dimensional case based on first
% order discretization. 
clear
clc
n=3;
M=20;
T=3;
points=multidimension(n,M);
P=[0.9 0.05 0.05; 0.05 0.9 0.05 ; 0.05 0.05 0.9];
C=[5 5 5 1];
obs=size(P,2);
point_num=size(points,2);
trans=zeros(point_num,point_num);
PTS=zeros(n+1,n);
Lambda=zeros(n+1,1);
V=zeros(T,point_num);
VC=zeros(T,point_num);
G=zeros(T,point_num);

p_temp=zeros(n,1);
dist=zeros(point_num,1);
for i=1:point_num
    for j=1:obs
        for k=1:n
            p_temp(k)=P(k,j)*points(k,i)/(P(:,j)'*points(:,i)); % using bayes rule to get p+
        end
        for ii=1:point_num
            dist(ii)=mydistance(p_temp,points(:,ii)); % find the distance between p+ and all grids 
        end
        [value, index]=sort(dist); % find the n nearest ones
        for jj=1:n  % p_temp=lambda_1*point1+lambda_2*point2+¡­¡­lambda*pointn ==> p_temp=[point1 point2 ... pointn]*[lambda_1; lambda_2;...lambda_n]
            PTS(:,jj)=[points(:,index(jj));1];
        end
        Lambda=pinv(PTS)*[p_temp;1];
        for kk=1:n
            trans(i,index(kk))=trans(i,index(kk))+Lambda(kk)*P(:,j)'*points(:,i); % update transition probability
        end
%         a=index(1);
%         b=index(2);
%         c=index(3);
%         trans(i,a)=trans(i,a)+lambda_a*P(:,j)'*points(:,i); % update transition probability
%         trans(i,b)=trans(i,b)+lambda_b*P(:,j)'*points(:,i);
%         trans(i,c)=trans(i,c)+lambda_c*P(:,j)'*points(:,i);
    end
end

% Value function at terminal step
Vn=zeros(1,n);
for i=1:point_num
    for j=1:n
        Vn(j)=points(j,i)*C(j);
    end
    [V(T,i),G(T,i)]=min(Vn);
end

% Value function
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


for i=1:point_num
    if G(1,i)==1
       plot3(points(1,i),points(2,i),points(3,i),'ro');
    else if G(1,i)==2
           plot3(points(1,i),points(2,i),points(3,i),'bo'); 
           hold on
        else if G(1,i)==3
                plot3(points(1,i),points(2,i),points(3,i),'ko');
                hold on
            else
                plot3(points(1,i),points(2,i),points(3,i),'go');
                hold on
            end
        end
    end
end