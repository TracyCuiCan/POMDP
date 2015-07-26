% This function finds the decentralized optimal thresholds using orthogonal
% search when the probability of observations is normally distributed.
function bnd = iterations_norm()
iter=100;
%Bound1=zeros(iter+1,2);
Bound1(1,:)=[0.1,0.9];
%Bound2=zeros(iter+1,2);
N = 10^3;
% P1 = [0.25 0.75; 0.6 0.4];
% P2 = [0.3 0.7; 0.8 0.2];
% C = [20, 60, 1]; % L1, L2, c
% P1 = [0.25 0.75; 0.6 0.4];
% P2 = [0.85 0.15; 0.35 0.65];
% C = [20 40 1];

x = linspace(-2, 2, 100);
f00 = normpdf(x, 0, 1);
f01 = normpdf(x, 1, 1);
f10 = normpdf(x, 0, 2);
f11 = normpdf(x, 1, 2);
P1 = [f00/sum(f00); f01/sum(f01)];
P2 = [f10/sum(f10); f11/sum(f11)];
C = [20 50 1];

for i=2:iter
    % set the threshold for sensor 1, optimize over threshold 2
    P1_u = absorption(P1,Bound1(i - 1,:));
    Bound2(i,:) = find_bound(P2, P1_u, C);
    % Use the result from last iteration, optimize over threshold 1
    P2_u = absorption(P2,Bound2(i,:));
    Bound1(i,:) = find_bound(P1, P2_u, C);
    
    if Bound2(i,1) == Bound2(i - 1,1) && Bound2(i,2) == Bound2(i - 1,2)% if bound remain the same, stop iteration
        break
    end
    
end
bnd = [Bound1(end,:), Bound2(end,:)];
end