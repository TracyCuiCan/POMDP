% This shows the discretized grids for three dimension. 
clear
clc
n=3;
M=20;
for i=1:n
    A(i,i)=1;
end
for i=1:n-1
    A(i,i+1)=-1;
end
B=(1/M)*A;

G_prime=zeros(n,(M+2)*(M+1)/2);
G=zeros(n,(M+2)*(M+1)/2);

G_prime(1,:)=M;
j=1;
for i=M:-1:0
    G_prime(2,j:(j+i))=i;
    for k=0:i
        G_prime(3,j+k)=k;
    end
    j=j+i+1;
end

G=B*G_prime;
x=G(1,:);
y=G(2,:);
z=G(3,:);
plot3(x,y,z,'o');
grid on;
p = {'linestyle','none','marker','o'};
line(x,y,z,p{:});

