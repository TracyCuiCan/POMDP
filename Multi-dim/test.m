%% This is to test whether multi-dimension function is working correctly
%%
%first_order_grid
%%
clear 
[V,G]=Discretize0(2, 1000, 3, [0.25 0.75; 0.6 0.4],[20 20 1]);
plot(0:0.001:1,fliplr(V(1,:)),'g')
%%
clear
[V,G]=Discretize1(2, 1000, 2, [0.25 0.75; 0.6 0.4],[20 20 1]);
plot(0:0.001:1,fliplr(V(1,:)),'g')