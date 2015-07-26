function bnd = iterations()
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
global P1 P2 C

for i=2:iter
    % set the threshold for sensor 1, optimize over threshold 2
    P1_u = absorption(P1,Bound1(i - 1,:));
%     for j = 1:N+1
%         P_u(j,:) = MonteCarlo(points(j), P1, Bound);
%     end
    Bound2(i,:) = find_bound(P2, P1_u, C);
    % Use the result from last iteration, optimize over threshold 1
    P2_u = absorption(P2,Bound2(i,:));
%     for j = 1:N+1
%         P_u(j,:) = MonteCarlo(points(j), P2, Bound);
%     end
    Bound1(i,:) = find_bound(P1, P2_u, C);
    
    if Bound2(i,1) == Bound2(i - 1,1) && Bound2(i,2) == Bound2(i - 1,2)% if bound remain the same, stop iteration
        break
    end
    
end
bnd = [Bound1(end,:), Bound2(end,:)];
end