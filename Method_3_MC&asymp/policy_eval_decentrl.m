% policy evaluation by sampling the Markov process of decentralized hypothesis testing
function value=valueat_des(bound)
P1 = [0.25 0.75; 0.6 0.4]; %[Y|h0; Y|h1]
P2 = [0.75 0.25; 0.4 0.6];
S=100;  % Times of sampling
N=1000;
cc=0.1;%cost to continue
L1=20; % the cost of making one error decision
L2=40; % the cost of making two error decisions
V=zeros(1,N+1);
result=zeros(size(P1,1),S);
points=linspace(0,1,N+1);

lowerbound1 = bound(1);
upperbound1 = bound(2);
lowerbound2 = bound(3);
upperbound2 = bound(4);

if (lowerbound1 < upperbound1)&&(lowerbound2 < upperbound2)
    for n=2:N+1
        p=points(n);% a prior probability of hypothesis H=0
        
        %sampling for sensor 1
        
        obs_pool=[randsrc(1,S*50,[1:size(P1,2);P1(1,:)]);
                   randsrc(1,S*50,[1:size(P1,2);P1(2,:)])];
        for h=1:size(P1,1) % p(u=0|H=0) and p(u=0|H=1)
            for s=1:S
                pi=p;
                count = 0;
                for i=((s-1)*50+1):s*50
                    index=obs_pool(h,i);
                    pi=pi*P1(1,index)/(pi*P1(1,index)+(1-pi)*P1(2,index)); % updating pi
                    count=count +1;
                    if pi<lowerbound1
                        result(h,s)=1;
                        break
                    else if pi>upperbound1
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
        % average length of observations required to reach a decision given H0 H1
        lgth=mean(obs_length,2);
        lgth11 = lgth(1);
        lgth12 = lgth(2);
        error = mean(result,2); % error(1)=p(u1|h0) error(2)=p(u1|h1)
        error11 = error(1);
        error12 = 1 - error(2);
        
        %sampling for sensor 2
        
        obs_pool=[randsrc(1,S*50,[1:size(P2,2);P2(1,:)]);
                   randsrc(1,S*50,[1:size(P2,2);P2(2,:)])];
        for h=1:size(P2,1) % p(u=0|H=0) and p(u=0|H=1)
            for s=1:S
                pi=p;
                count = 0;
                for i=((s-1)*50+1):s*50
                    index=obs_pool(h,i);
                    pi=pi*P2(1,index)/(pi*P2(1,index)+(1-pi)*P2(2,index)); % updating pi
                    count=count +1;
                    if pi<lowerbound2
                        result(h,s)=1;
                        break
                    else if pi>upperbound2
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
        lgth=mean(obs_length,2);
        lgth21 = lgth(1);
        lgth22 = lgth(2);
        error = mean(result,2); % error(1)=p(u1|h0) error(2)=p(u1|h1)
        error21 = error(1);
        error22 = 1 - error(2);
        
        V(n-1)=cc*(p*(lgth11+lgth21)+(1-p)*(lgth12+lgth22))+(1-p)*L1*(error12*(1-error22)+(1-error12)*error22)+p*L1*(error11*(1-error21)+(1-error11)*error21)+(1-p)*L2*error22*error12+p*L2*error11*error21;% We want to find the maximum value function using fminsearch
    end
    %plot(points,V)
    value=mean(V);
else
    value=inf;
end


end