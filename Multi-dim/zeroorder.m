% This gives optimal bounds for three dimensional case based on zeroth
% order discretization. 
clear
clc
n=3;
M=20;
T=4;
points=multidimension(n,M);
P=[0.9 0.05 0.05; 0.05 0.9 0.05 ; 0.05 0.05 0.9];
L=[5 5 5];
c=1;
obs=size(P,2);
point_num=size(points,2);
trans=zeros(point_num,point_num);
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
            dist(ii)=mydistance(p_temp,points(:,ii)); % find the closest point 
        end
        [dis, k]=min(dist);
        trans(i,k)=trans(i,k)+P(:,j)'*points(:,i); % update transition probability
    end
end

for i=1:point_num
    V1=points(1,i)*L(1);
    V2=points(2,i)*L(2);
    V3=points(3,i)*L(3);
    [V(T,i),G(T,i)]=min([V1 V2 V3]);
end

for t=T-1:-1:1
    for i=1:point_num
        Q=trans(i,:)*V(t+1,:)';
        VC(t,i)=c+Q;
        V1=points(1,i)*L(1);
        V2=points(2,i)*L(2);
        V3=points(3,i)*L(3);
        [V(t,i),G(t,i)]=min([V1 V2 V3 VC(t,i)]);
    end
end

figure
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
        
        
