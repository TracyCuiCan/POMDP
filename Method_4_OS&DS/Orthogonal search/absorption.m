function p_u = absorption(P,bound)
N = 10^3;
points = linspace(0,1,N+1);
lower = bound(1)*N+1 ; % index
upper = bound(2)*N+1 ;
interval = upper - lower + 1;

% transition probability
trans0=zeros(N+1,N+1);
trans1=zeros(N+1,N+1);
for i=1:length(points)
    beta=points(i);
        for j=1:size(P,2)
            p_temp=beta*P(1,j)/(beta*P(1,j)+(1-beta)*P(2,j)); % updating p
            k=max(1,ceil(p_temp*N)); % check which part the posterior lies in
            lambda=k-N*p_temp;
            % P(i,j|H = 0)
            p=P(1,j)*beta;
            trans0(i,k)=trans0(i,k)+p*lambda;
            trans0(i,k+1)=trans0(i,k+1)+p*(1-lambda);
            % P(i,j|H=1)
            p=P(2,j)*(1-beta); 
            trans1(i,k)=trans1(i,k)+p*lambda;
            trans1(i,k+1)=trans1(i,k+1)+p*(1-lambda);
        end
    trans0(i,:)=trans0(i,:)/sum(trans0(i,:));
    trans1(i,:)=trans1(i,:)/sum(trans1(i,:));
end
trans0(1,:) = 0; trans0(1,1) = 1;
trans1(N+1,:) = 0; trans1(N+1,N+1) = 1;


R = sum(trans0(lower:upper,1:lower - 1),2);% P(u = 1 | H0, bound)

tmp = trans0(lower:upper,lower:upper);

INV = eye(interval) - tmp; % N matrix for probability conditioned on H0
B = INV\R;
D1 = zeros(N+1,1);
D1(1:lower-1,1) = 1;
D1(lower:upper) = B;

% P(u = 0 | H1, bound)
R = sum(trans1(lower:upper,upper+1:N+1),2);
tmp = trans1(lower:upper,lower:upper);

INV = eye(interval) - tmp; % N matrix for probability conditioned on H1
B = INV\R;
D0 = zeros(N+1,1);
D0(upper+1:N+1,1) = 1;
D0(lower:upper) = B;

p_u = [1-D1, D1,D0,1-D0]; %u0|h0, u1|h0, u0|h1, u1|h1


