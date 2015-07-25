% policy evaluation by sampling the Markov process of centralized hypothesis testing
function value=valueat(bound)
P=[0.25 0.75; 0.6 0.4];
S=1000;  % Number of samplings
N=1000;
cc=0.01;%cost to continue
L0=20; % the cost of testing H0 when H1 is true
L1=20; % the cost of testing H1 when H0 is true
V=zeros(1,N+1);
result=zeros(size(P,1),S);
points=linspace(0,1,N+1);

lowerbound=bound(1);
upperbound=bound(2);

if (lowerbound<upperbound)
    for n=1:N+1
        p=points(n);% a prior probability of hypothesis H=0
        obs_pool=[randsrc(1,S*50,[1:size(P,2);P(1,:)]);
                   randsrc(1,S*50,[1:size(P,2);P(2,:)])];
        for h=1:size(P,1) % p(u=0|H=0) and p(u=0|H=1)
            %obs_pool=randsrc(1,S*50,[1:size(P,2);P(h,:)]);
            for s=1:S
                pi=p;
                count = 0;
                for i=((s-1)*50+1):s*50
                    index=obs_pool(h,i);
                    pi=pi*P(1,index)/(pi*P(1,index)+(1-pi)*P(2,index)); % updating pi
                    count=count +1;
                    if pi<lowerbound
                        result(h,s)=1;
                        break
                    else if pi>upperbound
                            result(h,s)=0;
                            break
                        else
                            continue
                        end
                    end
                end
                obs_length(h,s) = count;
            end
        end
        lgth=mean(obs_length,2); % average length of observations required to reach a decision given H0 H1
        prob=1-mean(result,2); % prob(1)=p(u0|h0) prob(2)=p(u0|h1)
        V(n)=p*(cc*lgth(1)+L0*(1-prob(1)))+(1-p)*(cc*lgth(2)+L1*prob(2));% We want to find the maximum value function using fminsearch
    end
    value=mean(V);
else
    value=inf;
end


end

