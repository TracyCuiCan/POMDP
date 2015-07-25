% This function is to use Monte Carlo Simulation to solve
% P(u=0|H=0,bounds), because the complete solution is computationally
% difficult to find due to the curse of dimension. Inputs are lower bound,
% upper bound and probability matrix, outputs are P(u|H,bounds)

% P=[0.25 0.75; 0.6 0.4]; % Probility matrix P(Y=y_i|H_i); row 1 H_0;
% clmn 1 Y_1
% tau=[0.3 0.8]; % lower and upper bound 


function [obs1,obs2, error1,error2]=MonteCarlo2(p,P1,P2,bound1,bound2)
tau1=bound1;
tau2=bound2;
N=10;  % Times of sampling
result1=zeros(2,N); % for two sensors and each sensor has two hypothesis
result2=zeros(2,N);
obs_length1=zeros(2,N);
obs_length2=zeros(2,N);


for h=1:2 % p(u=0|H=0) and p(u=0|H=1)
   for n=1:N
        obs_seq=randsrc(1,100,[1:size(P1,2);P1(h,:)]);% generate random sequence of observation index according to probability
        %sampling for sensor 1
        pi=p;
        for i=1:100
            index=obs_seq(i);
            pi=pi*P1(1,index)/(pi*P1(1,index)+(1-pi)*P1(2,index)); % updating pi
            if pi<tau1(1)
                result1(h,n)=1;
                break
            else if pi>tau1(2)
                    result1(h,n)=0;
                    break
                else
                    continue
                end
            end
        end
        obs_length1(h,n)=i;
        %sampling for sensor 2
        obs_seq=randsrc(1,100,[1:size(P2,2);P2(h,:)]);
        pi=p;
        for i=1:100
            index=obs_seq(i);
            pi=pi*P2(1,index)/(pi*P2(1,index)+(1-pi)*P2(2,index)); % updating pi
            if pi<tau2(1)
                result2(h,n)=1;
                break
            else if pi>tau2(2)
                    result2(h,n)=0;
                    break
                else
                    continue
                end
            end
        end
        obs_length2(h,n)=i;
    end
end
%sampling result for sensor 1
obs1=mean(obs_length1,2); % average length of observations required to reach a decision given H0 H1
error1=1-mean(result1,2); % by taking mean, we are actually calculating p(u=1|H)
conf1=[var(result1(1,:)); var(result1(2,:))]/sqrt(N); % by taking variance, we know our confidence in the simulation
%sampling result for sensor 2
obs2=mean(obs_length2,2); % average length of observations required to reach a decision given H0 H1
error2=1-mean(result2,2); % by taking mean, we are actually calculating p(u=1|H)
conf2=[var(result2(1,:)); var(result2(2,:))]/sqrt(N); % by taking variance, we know our confidence in the simulation

