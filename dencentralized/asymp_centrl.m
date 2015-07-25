% policy evaluation using asymptotic results for centralized setup
function value=asymp_centrl(bound)
% P1=[0.25 0.75; 0.6 0.4];
% P2 = [0.75 0.25; 0.4 0.6];
N=1000;
cc=0.01;%cost to continue
L0=20; % the cost of testing H0 when H1 is true
L1=20; % the cost of testing H1 when H0 is true
V=zeros(1,N+1);
points=linspace(0,1,N+1);
lowerbound=bound(1);
upperbound=bound(2);


if (lowerbound<upperbound)
    for n=2:N
       p=points(n);
       error1=lowerbound*(1-p)/p/(1-lowerbound); % P = (u1 | h0)
       error2=(1-upperbound)*p/(1-p)/upperbound; % P = (u0 | h1)
       lgth1=log(error2)/(log(0.6/0.25)*0.25 + log(0.4/0.75)*0.75);  %E(N|h0) probability matrix P1
       lgth2=log(1/error1)/ (log(0.6/0.25)*0.6 + log(0.4/0.75)*0.4); %E(N|h1) probability matrix P1
        
       V(n)=p*(cc*lgth1+L0*error1)+(1-p)*(cc*lgth2+L1*error2);% We want to find the maximum value function using fminsearch
    end
    value=mean(V);
else
    value=10^6;
end

