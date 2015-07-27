% This function is to generate a permutation {M=q1>=q2>=q3...}, n is the dimension, M is the desired parts

function G=multidimension(n,M)
A=npermutek([0:M],n); % call function npermutek to generate every possible permutation
%% delete those permutation no with descending order
for i=1:n-1
    temp1(:,i)=A(:,i)-A(:,i+1);
end
id=sum(temp1<0,2);
A(id>0,:)=[];

%% delete those permutation that doesn't start with M
id=find(A(:,1)~=M);
A(id,:)=[];

G_prime=A'; % each column of G_prime is one permutation that satisfies {M=q1>=q2>=q3...}

%% generate transition matrix B 
temp2=eye(n);
for i=1:n-1
    temp2(i,i+1)=-1; 
end
B=(1/M)*temp2;

%% G is the discretized points
G=B*G_prime;
end



        
        