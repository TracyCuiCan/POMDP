% This function is to compare the performance of alpha-vector, zeroth-order
% discretization and first-order discretization. 
load avec.mat;
load gridszero.mat;
load gridsfirst.mat;

figure
x=linspace(0,1,1000);
plot(x,gridszero,'b');
hold on;
y=linspace(0,1,1001);
plot(y,gridsfirst,'r');
hold on;
for i=1:size(avec,1)
    plot([0,1],[avec(i,1),avec(i,2)],'g');
    axis([0 1 0 7.5]);
    hold on;
end
xlabel('Belief Space');
ylabel('Value Function');
legend('0th order','1st order','\alpha-vector');