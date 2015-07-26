function [answer, fval] = fmin2()
options=optimset('MaxFunEvals',inf,'MaxIter',inf,'TolX',10^-3,'TolFun',10^-3);
N = 10^3;
global trans1 trans10 trans11 trans2 trans20 trans21 P1 P2 points
% P1 = [0.25 0.75; 0.6 0.4];
% P2 = [0.85 0.15; 0.35 0.65];

%points = linspace(0,1,N+1);

%% transition probability for DM1
trans1=zeros(N+1,N+1);
trans10=zeros(N+1,N+1);
trans11=zeros(N+1,N+1);
for i=1:length(points)
    beta=points(i);
        for j=1:size(P1,2)
            p_temp=beta*P1(1,j)/(beta*P1(1,j)+(1-beta)*P1(2,j)); % updating p
            k=max(1,ceil(p_temp*N)); % check which part the posterior lies in
            lambda=k-N*p_temp;
            % P(i,j)
            p=P1(1,j)*beta+P1(2,j)*(1-beta); % linear interpolation
            trans1(i,k)=trans1(i,k)+p*lambda;
            trans1(i,k+1)=trans1(i,k+1)+p*(1-lambda);
            % P(i,j|H = 0)
            p=P1(1,j);
            trans10(i,k)=trans10(i,k)+p*lambda;
            trans10(i,k+1)=trans10(i,k+1)+p*(1-lambda);
            % P(i,j|H=1)
            p=P1(2,j); 
            trans11(i,k)=trans11(i,k)+p*lambda;
            trans11(i,k+1)=trans11(i,k+1)+p*(1-lambda);
        end
    trans1(i,:)=trans1(i,:)/sum(trans1(i,:));
    trans10(i,:)=trans10(i,:)/sum(trans10(i,:));
    trans11(i,:)=trans11(i,:)/sum(trans11(i,:));
end
trans10(1,:) = 0; trans10(1,1) = 1;
trans11(N+1,:) = 0; trans11(N+1,N+1) = 1;
trans1 = sparse(trans1);
trans10 = sparse(trans10);
trans11 = sparse(trans11);

%% transition probability for DM2
trans2=zeros(N+1,N+1);
trans20=zeros(N+1,N+1);
trans21=zeros(N+1,N+1);
for i=1:length(points)
    beta=points(i);
        for j=1:size(P2,2)
            p_temp=beta*P2(1,j)/(beta*P2(1,j)+(1-beta)*P2(2,j)); % updating p
            k=max(1,ceil(p_temp*N)); % check which part the posterior lies in
            lambda=k-N*p_temp;
            % P(i,j)
            p=P2(1,j)*beta+P2(2,j)*(1-beta); % linear interpolation
            trans2(i,k)=trans2(i,k)+p*lambda;
            trans2(i,k+1)=trans2(i,k+1)+p*(1-lambda);
            % P(i,j|H = 0)
            p=P2(1,j);
            trans20(i,k)=trans20(i,k)+p*lambda;
            trans20(i,k+1)=trans20(i,k+1)+p*(1-lambda);
            % P(i,j|H=1)
            p=P2(2,j); 
            trans21(i,k)=trans21(i,k)+p*lambda;
            trans21(i,k+1)=trans21(i,k+1)+p*(1-lambda);
        end
    trans2(i,:)=trans2(i,:)/sum(trans2(i,:));
    trans20(i,:)=trans20(i,:)/sum(trans20(i,:));
    trans21(i,:)=trans21(i,:)/sum(trans21(i,:));
end
trans20(1,:) = 0; trans20(1,1) = 1;
trans21(N+1,:) = 0; trans21(N+1,N+1) = 1;
trans2 = sparse(trans2);
trans20 = sparse(trans20);
trans21 = sparse(trans21);

%% 
[answer, fval]=fminsearchbnd(@absorption2,[0.5 0.5 0.5 0.5],[0.001 0.001 0.001 0.001],[0.999 0.999 0.999 0.999],options);
end