% This function is doing policy evaluation that returns the performance of
% a specific bound for single decision maker.

function value = absorption1(bound)
N = 10^3;
c = 1;
L10 = 20;
L01 = 20;
%bound = [0.305 0.705];
points = linspace(0,1,N+1);
lower = round(bound(1),3)*N+1 ; % index
upper = round(bound(2),3)*N+1 ;
interval = upper - lower + 1;
global trans trans0 trans1 % These transition matrices are computed once in fmin.m to speed the computation
% computing E(N|bound)

tmp = trans(lower:upper,lower:upper);

lgth = zeros(1,N+1);
INV = eye(interval) - tmp; % N matrix for whole probability
lgth(lower:upper) = INV\ones(interval,1);

R = sum(trans0(lower:upper,1:lower - 1),2);% P(u = 1 | H0, bound)

tmp = trans0(lower:upper,lower:upper);

INV = eye(interval) - tmp; % N matrix for probability conditioned on H0
B = INV\R;
D1 = zeros(1,N+1);
D1(1,1:lower-1) = 1;
D1(lower:upper) = B;

% P(u = 0 | H1, bound)
R = sum(trans1(lower:upper,upper+1:N+1),2);
tmp = trans1(lower:upper,lower:upper);

INV = eye(interval) - tmp; % N matrix for probability conditioned on H1
B = INV\R;
D0 = zeros(1,N+1);
D0(1,upper+1:N+1) = 1;
D0(lower:upper) = B;

%computing value function
for i = 1:N+1
    V(i) = lgth(i)*c + D0(i)*(1-points(i))*L01 + D1(i)*points(i)*L10;
end

value = mean(V);



