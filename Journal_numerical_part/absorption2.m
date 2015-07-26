function value = absorption2(bound)
global trans1 trans10 trans11 trans2 trans20 trans21 C points
%C = [20 50 1];
L1 = C(1);
L2 = C(2);
c = C(3);
N = 10^3;
%points = linspace(0,1,N+1);

%% DM1
lower1 = round(bound(1),3)*N+1 ; % index
upper1 = round(bound(2),3)*N+1 ;
inter1 = upper1 - lower1 + 1;
% computing E(N|bound)

tmp = trans1(lower1:upper1,lower1:upper1);

lgth1 = zeros(1,N+1);
INV = eye(inter1) - tmp; % N matrix for whole probability
lgth1(lower1:upper1) = INV\ones(inter1,1);

R = sum(trans10(lower1:upper1,1:lower1 - 1),2);% P(u = 1 | H0, bound)

tmp = trans10(lower1:upper1,lower1:upper1);

INV = eye(inter1) - tmp; % N matrix for probability conditioned on H0
B = INV\R;
D11 = zeros(1,N+1);
D11(1,1:lower1-1) = 1;
D11(lower1:upper1) = B;

R = sum(trans11(lower1:upper1,upper1+1:N+1),2); % P(u = 0 | H1, bound)
tmp = trans11(lower1:upper1,lower1:upper1);

INV = eye(inter1) - tmp; % N matrix for probability conditioned on H1
B = INV\R;
D10 = zeros(1,N+1);
D10(1,upper1+1:N+1) = 1;
D10(lower1:upper1) = B;

%% DM2
lower2 = round(bound(3),3)*N+1 ; % index
upper2 = round(bound(4),3)*N+1 ;
inter2 = upper2 - lower2 + 1;
% computing E(N|bound)

tmp = trans2(lower2:upper2,lower2:upper2);

lgth2 = zeros(1,N+1);
INV = eye(inter2) - tmp; % N matrix for whole probability
lgth2(lower2:upper2) = INV\ones(inter2,1);

R = sum(trans20(lower2:upper2,1:lower2 - 1),2);% P(u = 1 | H0, bound)

tmp = trans20(lower2:upper2,lower2:upper2);

INV = eye(inter2) - tmp; % N matrix for probability conditioned on H0
B = INV\R;
D21 = zeros(1,N+1);
D21(1,1:lower2-1) = 1;
D21(lower2:upper2) = B;

R = sum(trans21(lower2:upper2,upper2+1:N+1),2); % P(u = 0 | H1, bound)
tmp = trans21(lower2:upper2,lower2:upper2);

INV = eye(inter2) - tmp; % N matrix for probability conditioned on H1
B = INV\R;
D20 = zeros(1,N+1);
D20(1,upper2+1:N+1) = 1;
D20(lower2:upper2) = B;

%% computing value function
for i = 1:N+1
    V(i) = (lgth1(i)+lgth2(i))*c + (D11(i)*D21(i)*points(i)+D10(i)*D20(i)*(1-points(i)))*L2 + ((D11(i)*(1-D21(i))+D21(i)*(1-D11(i)))*points(i)+(D10(i)*(1-D20(i))+D20(i)*(1-D10(i)))*(1-points(i)))*L1;
end

value = mean(V,2);

