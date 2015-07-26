iter = 100;
val = zeros(iter,2); % first column OS, second column DS
thrshd = cell(iter,2); % first column OS, second column DS
N = 10^3;
ds_time = zeros(1,iter);
os_time = zeros(1,iter);
%tic
for i = 1:iter
    global P1 P2 C points
    points=linspace(0,1,N+1);
    %points = [0,unique(rand(1,N-1)),1];
%         P1 = [0.25 0.75; 0.6 0.4];
%         P2 = [0.8 0.2; 0.3 0.7];
%         C = [20 50 1];
    
    tmp1 = rand(2,2);
    P1 = [tmp1(1,:)/sum(tmp1(1,:)); tmp1(2,:)/sum(tmp1(2,:))];
    tmp1 = rand(2,2);
    P2 = [tmp1(1,:)/sum(tmp1(1,:)); tmp1(2,:)/sum(tmp1(2,:))];
    C=[20,50,1]; % the cost matrix 
    
    %C=[randi([2,40]),randi([40,80]),1];
    %     temp = rand(2);
    %     P1 = [temp(1) 1- temp(1); temp(2) 1- temp(2)]; % Probability matrix for two hypothesis
    %     temp = rand(2);
    %     P2 = [temp(1) 1- temp(2,1); temp(2) 1- temp(2)];
    t1 = clock;
    [thrshd{i, 2}, val(i,2)] = fmin2();
    t2 = clock;
    ds_time(i) = etime(t2,t1);
    
    t3 = clock;
    thrshd{i,1} = iterations();
    val(i,1) = absorption2(thrshd{i,1});
    t4 = clock;
    os_time(i) = etime(t4,t3);
    %[thrshd{i,3}, val(i,3)] = first_grid();
    clear var tmp1 tmp2
    clear global var P1 P2 C trans1 trans10 trans11 trans2 trans20 trans21
end
%toc