% policy evaluation using asymptotic results for decentralized setup
function value=asymp_decentrl(bound)
%bound =[bound_sensor1, bound_sensor_2];
% P1=[0.25 0.75; 0.6 0.4];
% P2=[0.75 0.25; 0.4 0.6];
N=1000;
cc=0.01;%cost to continue
L1=10; % the cost of making one wrong decision 
L2=20; % the cost of making two wrong decisions
V=zeros(1,N-1);
points=linspace(0,1,N+1);
LB1=bound(1);
UB1=bound(2);
LB2=bound(3);
UB2=bound(4);


if (LB1>UB1)&&(LB2>UB2)
    value=inf;
else
    for n=2:N
        p=points(n);
        %error and observations for sensor 1
       error11=LB1*(1-p)/p/(1-LB1); %p(u=1|h0)
       error12=(1-UB1)*p/(1-p)/UB1; %p(u=0|h1)
       lgth11=log(error12)/(log(0.6/0.25)*0.25 + log(0.4/0.75)*0.75);  %E(N|h0) probability matrix P1
       lgth12=log(1/error11)/ (log(0.6/0.25)*0.6 + log(0.4/0.75)*0.4); %E(N|h1) probability matrix P1
%        V(n-1)=p*(cc*lgth11+L0*error11)+(1-p)*(cc*lgth12+L1*error12);
       %error and observations for sensor 2
       error21=LB2*(1-p)/p/(1-LB2);
       error22=(1-UB2)*p/(1-p)/UB2;
       lgth21=log(error22)/(log(0.6/0.25)*0.25 + log(0.4/0.75)*0.75);  %E(N|h0) probability matrix P2
       lgth22=log(1/error21)/ (log(0.6/0.25)*0.6 + log(0.4/0.75)*0.4); %E(N|h1) probability matrix P2
%       V(n-1)=p*(cc*lgth21+L0*error21)+(1-p)*(cc*lgth22+L1*error22);
       V(n-1)=cc*(p*(lgth11+lgth21)+(1-p)*(lgth12+lgth22))+(1-p)*L1*(error12*(1-error22)+(1-error12)*error22)+p*L1*(error11*(1-error21)+(1-error11)*error21)+(1-p)*L2*error22*error12+p*L2*error11*error21;% We want to find the maximum value function using fminsearch
    end
    value=mean(V);
end

end