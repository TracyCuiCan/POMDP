% This is the test case for randomly chosen probability matrix P1 and P2. 
tic
t = 1;
val = zeros(t,2); % first column OS, second column DS
thrshd = cell(t,2); % first column OS, second column DS
for i = 1:t
    global P1 P2 C
    
%     P1 = [0.25 0.75; 0.6 0.4];
%     P2 = [0.8 0.2; 0.3 0.7];
%     C = [20 50 1];

    tmp1 = rand(2,8);
    P1 = [tmp1(1,:)/sum(tmp1(1,:)); tmp1(2,:)/sum(tmp1(2,:))];
    tmp1 = rand(2,8);
    P2 = [tmp1(1,:)/sum(tmp1(1,:)); tmp1(2,:)/sum(tmp1(2,:))];
%     
%     C=[randi([2,40]),randi([40,80]),1]; % Generate the cost matrix
%     temp = rand(2);
%     P1 = [temp(1) 1- temp(1); temp(2) 1- temp(2)]; % Probability matrix for two hypothesis
%     temp = rand(2);
%     P2 = [temp(1) 1- temp(2,1); temp(2) 1- temp(2)];

    [thrshd{i, 2}, val(i,2)] = fmin2();
    thrshd{i,1} = iterations();
    val(i,1) = absorption2(thrshd{i,1});  
    %[thrshd{i,3}, val(i,3)] = first_grid();
    clear P1 P2 C temp
end
toc