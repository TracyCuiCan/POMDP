% This code is to testify whether function Sequential is working correctly.
% The return of test_Sequential is compared with the result of bounds using
% alpha-vectors;
clear
P=[0.25  0.75; 0.6 0.4];
C=[20 20 1];
horizon=3;
% The result should be taolower =[0.5000 0.391 0.326],taoupper =[0.5000
% 0.647 0.691]

[lower,upper]=Sequential(P,C,horizon);
lower;
upper;
