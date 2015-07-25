%This function is to generate alpha-vector for finite horizon scenario
% One can adapt this to different case by changing parameters.
clear;
clc;
tic;

P=[0.25 0.75; 0.6 0.4];
obs=size(P,1);
T=6; % horizon
cc=1; %cost to continue
L0=20; % the cost of testing H0 when H1 is true
L1=20; % the cost of testing H1 when H0 is true
x=linspace(0,1,10^3+1);

AVEC = cell(1,T);  % store all alpha-vectors
AVEC{1} = [0 L1; L0 0];
avecfirst = zeros(1, 10^3+1);
% when t=1, we only have first two term
for t=1
    Avec=[0 L1; L0 0]; 
    Number(t)=size(Avec,1);
    x0=intersec([0 L1],[1,0],[0 0],[1 L0]);
    y0=[0 L1]*[x0; 1-x0];
    figure(t)
    taolower(1)=x0;
    taoupper(1)=x0;
    for i=1:size(Avec,1)
        plot([0,1],[Avec(i,:)*[0;1],Avec(i,:)*[1;0]],'b','linewidth',2.5);
        axis([0 1 0 y0]);
        hold on;
        xx=['p(',num2str(t),')'];
        yy=['V(',num2str(t),')'];
        xlabel(xx);
        ylabel(yy);
        plot(taolower(t),L0*taolower(t),'r.','MarkerSize',15);
        t_text=['\alpha=\beta=',num2str(taolower(t))];
        text(taolower(1)-0.2,L0*taolower(t)-1,t_text);
    end
end

for t=2:T
    %find all alpha-vectors under different Ys
    for j=1:obs
        temp(:,1)=Avec(:,1)*P(1,j); % a*f0(y)
        temp(:,2)=Avec(:,2)*P(2,j); %b*f1(y)
        a(j,:)=reshape(temp',1,size(temp,1)*size(temp,2));
    end
    
    % find all possible summation of alpha-vectors under different
    % observation Ys
    
    A=(reshape(a(1,:),2,length(a(1,:))/2))';
    for j=2:obs
        B=(reshape(a(j,:),2,length(a(j,:))/2))';
        A=minksum(A,B);
    end
    C=unique(A,'rows');
    ind=~any(C,2);
    C(ind,:)=[];
    
    % find the cross points of those vectors
    d=node(C);
    d=unique(d);
    d=[0,sort(d),1];
    
    % find minimum corresponding alpha-vectors
    for i=1:length(d)-1
        dtemp=(d(i)+d(i+1))/2;
        Ctemp=C*[dtemp; 1-dtemp];
        [V,index]=min(Ctemp,[],1);
        avec(i,1)=C(index,:)*[0;1];
        avec(i,2)=C(index,:)*[1;0];
    end
    avec=unique(avec,'rows');
    ind=~any(avec,2);
    avec(ind,:)=[];
    avec=avec+cc;
    
    % find the upper and lower bound for each horizon
    n=size(avec,1);
    for i=1:n
        upper(i)=intersec([x0, y0],[1,0],[0,avec(i,1)],[1,avec(i,2)]);
        lower(i)=intersec([0, 0],[x0,y0],[0,avec(i,1)],[1,avec(i,2)]);
    end
    lower=lower(find(lower>0 & lower<1));
    upper=upper(find(upper>0 & upper<1));
    taolower(t)=min(lower);
    taoupper(t)=max(upper);
    
    %plot figures, T is the horizon
    figure(t)
    for i=1:size(avec,1)
        plot([0,1],[avec(i,1),avec(i,2)],'g','linewidth',1.5);
        axis([0 1 0 y0]);
        hold on;
    end
    plot([0,1],[L1,0],'g','linewidth',1.5);
    hold on
    plot([0,1],[0,L0],'g','linewidth',1.5)
    plot(taolower(t),20*taolower(t),'r.','MarkerSize',15);
    plot(taoupper(t),20*(1-taoupper(t)),'r.','MarkerSize',15);
    t_text_1=['\alpha=',num2str(taolower(t))];
    t_text_2=['\beta=',num2str(taoupper(t))];
    text(taolower(t)-0.1,20*taolower(t)-1,t_text_1);
    text(taoupper(t)-0.05,20*(1-taoupper(t))-1,t_text_2);
    xx=['p(',num2str(t),')'];
    yy=['V(',num2str(t),')'];
    xlabel(xx);
    ylabel(yy);
    hold on
    
    %transform Avec back into a row vector and clear storage
    Avec=[avec; 0 L1; L0 0];
    for i=1:length(x)
        value = zeros(size(Avec,1),1);
        for j = 1: size(Avec, 1)
            value(j) = (Avec(j,2) - Avec(j,1))*x(i)+Avec(j,1);
        end
        avecfirst(i) = min(value);
    end
    plot(x,avecfirst,'b.','MarkerSize',7);
    AVEC{t} = Avec;
    
    
    % Avec=[ 0 L1; L0 0;avec;];
    Number(t)=size(Avec,1);
    clear a d B C  avec temp
end

% figure(1)
% t=1:T;
% plot(t,taolower,'r');
% hold on
% plot(t,taoupper,'b');
% figure(2)
% plot(t,Number(t),'k');

taolower;
taoupper;
N = 10^3;
points=linspace(0,1,N+1);
for n = 1:N + 1
    for i = 1:size(Avec,1)
        V(n, i) = (Avec(i,2) - Avec(i,1))*points(n) + Avec(i,1);
    end
    W(n) = min(V(n,:));
end
value = mean(W);



toc

