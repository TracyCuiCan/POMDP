load('random500.mat');
os = (val(:,1) - val(:,2))./val(:,1); % var1
ds = (val(:,2) - val(:,1))./val(:,2); % var2
%[cnt, bins] = hist(os,100); % use cnt to check the counts

subplot(1,2,1);
hist(os,100);
xlabel('Relative Error');
ylabel('Counts');
title('Histogram of \Delta J_{OS}');
xlim([-1.1 0.6]); % put two histograms in same frame in order to compare
ylim([0 120]); 

subplot(1,2,2);
hist(ds,100);
xlabel('Relative Error');
ylabel('Counts');
title('Histogram of \Delta J_{DS}');
xlim([-1.1 0.6]); % put two histograms in same frame in order to compare
ylim([0 120]); 

