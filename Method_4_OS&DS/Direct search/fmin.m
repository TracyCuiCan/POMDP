% This function computes three transition matrices and finds the optimal
% threshold using direct search for centralized hypothesis testing. 

tic 
%options=optimset('Display','Iter','MaxFunEvals',inf,'MaxIter',inf);
options=optimset('MaxFunEvals',inf,'MaxIter',inf,'TolX',0.001,'TolX',0.001);
global trans trans0 trans1
N = 10^3;
points = linspace(0,1,N+1);
P=[0.25 0.75; 0.6 0.4];
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

% M = 100;
% count = 0;
% for i = 1:M-1
%     for j = (i+1):M
%          v = absorption([i/M j/M]);
%          plot3(i,j,v)
%          grid on
%          hold on
%     end
% end

[bound, fval, exitflag, output]=fminsearchbnd(@absorption1,[0.3 0.7],[0.001 0.001],[0.999 0.999],options)

toc